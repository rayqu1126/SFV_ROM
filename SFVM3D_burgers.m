%% SFV and ROM for 1D Burgers equations with 2D stochstic variables
clear
close all
addpath('func')

% Set default plot
set(0, 'defaultaxesfontsize',24,'defaultaxeslinewidth',2,...
       'defaultlinelinewidth',2,'defaultpatchlinewidth',2,...
       'defaulttextfontsize',24,'defaulttextinterpreter','latex');
set(0, 'DefaultFigurePosition',[0 0 600 400]);

% Initialize mesh (Nx: physical Ny: stochastic)
Nx = 64; 
Ny = 32;

xpts = linspace(0,1,Nx+1);
ypts = linspace(0,1,Ny+1);
zpts = linspace(0,1,Ny+1);

x = 0.5*(xpts(1:Nx)+xpts(2:Nx+1));
y = 0.5*(ypts(1:Ny)+ypts(2:Ny+1));
z = 0.5*(zpts(1:Ny)+zpts(2:Ny+1));

[Y,X,Z] = meshgrid(y,x,z);

% Initialize solution array
ics = sin(2*pi*X) + 0.5 * sin(2*pi*Y) + Z;

t0=0;
tmax=0.2;
tspan = [t0 tmax];

dx = 1/Nx;
dy = 1/Ny;
dz = 1/Ny;

% Assemble parameters
params.N = Ny;
params.Nx = Nx;
params.dx = dx;
params.dy = dy;
params.dz = dz;

% SFV - WENO with reconstructed states
options = odeset('RelTol',1e-6,'AbsTol',1e-8);

tic
[~,U_state] = ode45(@(t,U) rhs_3D_burgers_state(t,U,params), tspan, ics, options);
toc

sol_state =reshape(U_state(end,:),Nx,Ny,Ny);

% f1 = @() ode45(@(t,U) rhs_3D_burgers_state(t,U,params), tspan, ics, options);
% timeit(f1)

% SFV - WENO with reconstructed fluxes
tic
[~,U_flux] = ode45(@(t,U) rhs_3D_burgers_flux(t,U,params), tspan, ics, options);
toc 

sol_flux=reshape(U_flux(end,:),Nx,Ny,Ny);

% f2 = @() ode45(@(t,U) rhs_3D_burgers_flux(t,U,params), tspan, ics, options);
% timeit(f2)

% Compute relative difference
er_flux = sum(abs(sol_state - sol_flux)...
    * (dx*dy*dz),"all") / sum(abs(sol_state) * (dx*dy*dz),"all");

disp("The relative L1 difference (recosntructed states and reconstructed fluxes) is " + er_flux)


% Solution plots (with reconstructed fluxes)
[YY,XX] = meshgrid(y,x);
figure
surf(XX,YY,sol_flux(:,:,end))
xlabel('$x$'), ylabel('$y_1$'), zlabel('$U$')
view(52.5, 30)

figure
surf(XX,YY,reshape(sol_flux(:,end,:),Nx,Ny))
xlabel('$x$'), ylabel('$y_2$'), zlabel('$U$')
view(52.5, 30)

% Compute and plot mean and std of solution (with reconstructed fluxes)
sol_mean = sum(sum(sol_flux, 2), 3) * dy * dz;  
sol_diff = sol_flux - sol_mean;  
sol_var = sum(sum(sol_diff.^2, 2), 3) / (Ny^2 - 1); 
sol_std = sqrt(sol_var); 
% 
figure
h1 = plot(x, sol_mean, 'LineWidth', 2, 'Color', 'r'); 
hold on
h2 = plot(x, sol_mean + sol_std, 'LineWidth', 2, 'Color', 'r', 'LineStyle', '--');
plot(x, sol_mean - sol_std, 'LineWidth', 2, 'Color', 'r', 'LineStyle', '--');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$U$', 'Interpreter', 'latex');
legend([h1, h2], {'Mean', 'Mean $\pm$ Std'}, 'Interpreter', 'latex',   'Location', 'best');


%% ROM - POD
% 1D sampling
tspan_sample = linspace(t0, tmax, 200);
F_sample = zeros(4*Ny^2,Nx*length(tspan_sample));

% Coordinates y and z of quad points volume by volume
yquadL = y + 0.5*(-1/sqrt(3)) * dy;
yquadR = y + 0.5*(1/sqrt(3)) * dy;   
yquad = reshape([yquadL; yquadR; yquadL; yquadR], [], 1)';
yquad = repmat(yquad,1,Ny);

zquadL = z + 0.5*(-1/sqrt(3)) * dz;
zquadR = z + 0.5*(1/sqrt(3)) * dz;
zquad = repmat([zquadL; zquadL; zquadR; zquadR],Ny,1);
zquad = reshape(zquad,[],1)';

for i = 1:length(yquad)
    ics_sample = sin(2*pi*x') + 0.5 * sin(2*pi*yquad(i)) + zquad(i);
    F_sample(i,:) = sampling_burgers(Nx, ics_sample, tspan_sample);
end

% Construct reduced basis
Nmode = 10;
[Vn,s,~] = svd(F_sample);
Vn = Vn(:,1:Nmode);

% Integrals of bases
Bn = zeros(Ny^2,Nmode);
for i = 1:Ny^2
    Bn(i,:) = (Vn(4*i-3,:) +Vn(4*i-2,:)+Vn(4*i-1,:)+ Vn(4*i,:)) /4;
end

% Q-DEIM for hyper-reduction
[~,~,P] = qr(Vn',"vector");
NHR = Nmode;
ids = P(1:NHR);
Vn_HR = Vn(ids,:);

%
params.Vinv = pinv(Vn_HR);
params.B = Bn;
params.ids = ids;

tic
[~,U_ROM] = ode45(@(t,U) rhs_3D_burgers_interp(t,U,params), tspan, ics, options);
toc

sol_ROM = reshape(U_ROM(end,:),Nx,Ny,Ny);

% Compute ROM relative error
er_ROM = sum(abs(sol_ROM - sol_flux)...
    * (dx*dy*dz),"all") / sum(abs(sol_flux) * (dx*dy*dz),"all");

disp("The relative L1 error (reconstructed fluxes and ROM) is " + er_ROM)

% ROM solution plots
figure
surf(XX,YY,sol_ROM(:,:,end))
xlabel('$x$'), ylabel('$y_1$'), zlabel('$U$')
view(52.5, 30)

figure
surf(XX,YY,reshape(sol_ROM(:,end,:),Nx,Ny))
xlabel('$x$'), ylabel('$y_1$'), zlabel('$U$')
view(52.5, 30)
% 

% Compute and plot ROM mean and std
sol_mean = sum(sum(sol_ROM, 2), 3) * dy * dz;  
sol_diff = sol_ROM - sol_mean;  
sol_var = sum(sum(sol_diff.^2, 2), 3) / (Ny^2 - 1); 
sol_std = sqrt(sol_var); 
% 
figure
h1 = plot(x, sol_mean, 'LineWidth', 2, 'Color', 'r'); 
hold on
h2 = plot(x, sol_mean + sol_std, 'LineWidth', 2, 'Color', 'r', 'LineStyle', '--');
plot(x, sol_mean - sol_std, 'LineWidth', 2, 'Color', 'r', 'LineStyle', '--');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$U$', 'Interpreter', 'latex');
legend([h1, h2], {'Mean', 'Mean $\pm$ Std'}, 'Interpreter', 'latex',   'Location', 'best');
