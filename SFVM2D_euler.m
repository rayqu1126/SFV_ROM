%% SFVM and ROM for 1D Euler equations with 1D stochastic variables
clear
close all
addpath('func')

% Set default plot
set(0, 'defaultaxesfontsize',24,'defaultaxeslinewidth',2,...
       'defaultlinelinewidth',2,'defaultpatchlinewidth',2,...
       'defaulttextfontsize',24,'defaulttextinterpreter','latex');
set(0, 'DefaultFigurePosition',[0 0 600 400]);

% Initialize mesh (Nx: physical Ny: stochastic)
Nx = 128;
Ny = 32;

xpts = linspace(0,1,Nx+1);
ypts = linspace(0,1,Ny+1);

x = 0.5*(xpts(1:Nx)+xpts(2:Nx+1));
y = 0.5*(ypts(1:Ny)+ypts(2:Ny+1));

[Y,X] = meshgrid(y,x);

% Initialize solution array
t0=0;
tmax=0.2;
tspan = [t0 tmax];

u = zeros(Nx,Ny);
rho = ones(Nx,Ny) * 0.125;
p = ones(Nx,Ny) * 0.1;


rho(X < 0.475 + 0.05*Y) = 1; 
p(X < 0.475 + 0.05*Y) = 1;   
gamma = 1.4;
E = p / (gamma - 1) + 0.5 * rho .* u.^2;

ics = [rho; rho .* u; E];

% Assemble parameters
dx = 1/Nx;
dy = 1/Ny;

params.N = Ny;
params.Nx = Nx;
params.dx = dx;
params.dy = dy;
params.gamma = gamma;

% SFV - with reconstructed states
options = odeset('NonNegative', numel(ics), 'RelTol',1e-6,'AbsTol',1e-8);
tic
[~,U_state] = ode45(@(t,U) rhs_2D_euler_state(t,U,params), tspan, ics, options);
toc

sol_state = reshape(U_state(end,:),3*Nx,Ny);

% f1 = @() ode45(@(t,U) rhs_2D_euler_state(t,U,params), tspan, ics, options);
% timeit(f1)

% SFV - with reconstructed fluxes
tic
[~,U_flux] = ode45(@(t,U) rhs_2D_euler_flux(t,U,params), tspan, ics, options);
toc

% f2 = @() ode45(@(t,U) rhs_2D_euler_flux(t,U,params), tspan, ics, options);
% timeit(f2)

sol_flux = reshape(U_flux(end,:),3*Nx,Ny);
sol_flux_rho = sol_flux(1:Nx,:);
sol_flux_rhou = sol_flux(Nx+1:2*Nx,:);
sol_flux_u = sol_flux_rhou ./ sol_flux_rho;


% Compute relative difference
er_flux = sum(abs(sol_state - sol_flux) * (dx*dy),"all") ...
    / sum(abs(sol_state) * (dx*dy),"all");

disp("The relative L1 difference (recosntructed states and reconstructed fluxes) is " + er_flux)


% Solution plots (with reconstructed fluxes)
figure
surf(X,Y,sol_flux_rho)
xlabel('$x$'), ylabel('$y$'), zlabel('$\rho$')
zlim([-0.15 1.45])
view(52.5, 30)

figure
surf(X,Y,sol_flux_u)
xlabel('$x$'), ylabel('$y$'), zlabel('$u$')
zlim([-0.15 1.45])
view(52.5, 30)

% Compute and plot mean and std of solution (with reconstructed fluxes)
sol_rho_mean = mean(sol_flux_rho, 2);
sol_rho_diff = sol_flux_rho - sol_rho_mean;  
sol_rho_var = sum(sol_rho_diff.^2, 2) / (Ny - 1);  
sol_rho_std = sqrt(sol_rho_var);

figure
h1 = plot(x, sol_rho_mean, 'LineWidth', 2, 'Color', 'r'); 
hold on
h2 = plot(x, sol_rho_mean + sol_rho_std, 'LineWidth', 2, 'Color', 'r', 'LineStyle', '--');
plot(x, sol_rho_mean - sol_rho_std, 'LineWidth', 2, 'Color', 'r', 'LineStyle', '--');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$\rho$', 'Interpreter', 'latex');
legend([h1, h2], {'Mean', 'Mean $\pm$ Std'}, 'Interpreter', 'latex',   'Location', 'best');

%
sol_u_mean = mean(sol_flux_u, 2);
sol_u_diff = sol_flux_u - sol_u_mean;  
sol_u_var = sum(sol_u_diff.^2, 2) / (Ny - 1);  
sol_u_std = sqrt(sol_u_var);

figure
h1 = plot(x, sol_u_mean, 'LineWidth', 2, 'Color', 'r'); 
hold on
h2 = plot(x, sol_u_mean + sol_u_std, 'LineWidth', 2, 'Color', 'r', 'LineStyle', '--');
plot(x, sol_u_mean - sol_u_std, 'LineWidth', 2, 'Color', 'r', 'LineStyle', '--');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$u$', 'Interpreter', 'latex');
legend([h1, h2], {'Mean', 'Mean $\pm$ Std'}, 'Interpreter', 'latex',   'Location', 'NW');

%% ROM - POD
% 1D sampling
tspan_sample = linspace(t0, tmax, 100);
quadL = y + 0.5*(-1/sqrt(3)) * dy;
quadR = y + 0.5*(1/sqrt(3)) * dy;  

F_sample = zeros(2*Ny, 2*3*Nx*length(tspan_sample));
quad_combined = reshape([quadL; quadR], [], 1)';
for i = 1:length(quad_combined)
    u = zeros(Nx,1);
    rho = ones(Nx,1) * 0.125;
    p = ones(Nx,1) * 0.1;
    rho(x < 0.475+0.05*quad_combined(i)) = 1; 
    p(x < 0.475+0.05*quad_combined(i)) = 1;     
    E = p / (gamma - 1) + 0.5 * rho .* u.^2;    
    ics_sample = [rho; rho.*u; E];
    F_sample(i,:)= sampling_euler(Nx, gamma, ics_sample, tspan_sample);
end

% Construct reduced basis
Nmode = 40;
[Vn,s,~] = svd(F_sample,"econ"); % ECON-SVD for large matrix
Vn = Vn(:,1:Nmode);

% Integral of basis
Bn = zeros(Ny,Nmode);
for i = 1:Ny
    Bn(i,:) = (Vn(2*i-1,:) + Vn(2*i,:)) /2;
end

% Q-DEIM hyper-reduction
[~,~,P] = qr(Vn',"vector");

NHR = 2 * Ny;
ids = P(1:NHR);

Vn_HR = Vn(ids,:);

%
params.Vinv = pinv(Vn_HR);
params.B = Bn;
params.ids = ids;

tic
[~,U_ROM] = ode45(@(t,U) rhs_2D_euler_interp(t,U,params), tspan, ics, options);
toc
%

sol_ROM = reshape(U_ROM(end,:),3*Nx,Ny);
sol_ROM_rho = sol_ROM(1:Nx,:);
sol_ROM_rhou = sol_ROM(Nx+1:2*Nx,:);
sol_ROM_u = sol_ROM_rhou ./ sol_ROM_rho;

% Compute ROM relative error
er_ROM = sum(abs(sol_flux - sol_ROM) * (dx*dy),"all") ...
    / sum(abs(sol_flux) * (dx*dy),"all");

disp("The relative L1 error (reconstructed fluxes and ROM) is " + er_ROM)

% ROM solution plots
figure
surf(X,Y,sol_ROM_rho)
xlabel('$x$'), ylabel('$y$'), zlabel('$\rho$')
zlim([-0.15 1.45])
view(52.5, 30)

figure
surf(X,Y,sol_ROM_u)
xlabel('$x$'), ylabel('$y$'), zlabel('$u$')
zlim([-0.15 1.45])
view(52.5, 30)

% Compute and plot mean and std of ROM solutions
sol_rho_mean = mean(sol_ROM_rho, 2);
sol_rho_diff = sol_ROM_rho - sol_rho_mean;  
sol_rho_var = sum(sol_rho_diff.^2, 2) / (Ny - 1);  
sol_rho_std = sqrt(sol_rho_var);

figure
h1 = plot(x, sol_rho_mean, 'LineWidth', 2, 'Color', 'r'); 
hold on
h2 = plot(x, sol_rho_mean + sol_rho_std, 'LineWidth', 2, 'Color', 'r', 'LineStyle', '--');
plot(x, sol_rho_mean - sol_rho_std, 'LineWidth', 2, 'Color', 'r', 'LineStyle', '--');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$\rho$', 'Interpreter', 'latex');
legend([h1, h2], {'Mean', 'Mean $\pm$ Std'}, 'Interpreter', 'latex',   'Location', 'best');

% 
sol_u_mean = mean(sol_ROM_u, 2);
sol_u_diff = sol_ROM_u - sol_u_mean;  
sol_u_var = sum(sol_u_diff.^2, 2) / (Ny - 1);  
sol_u_std = sqrt(sol_u_var);

figure
h1 = plot(x, sol_u_mean, 'LineWidth', 2, 'Color', 'r'); 
hold on
h2 = plot(x, sol_u_mean + sol_u_std, 'LineWidth', 2, 'Color', 'r', 'LineStyle', '--');
plot(x, sol_u_mean - sol_u_std, 'LineWidth', 2, 'Color', 'r', 'LineStyle', '--');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$u$', 'Interpreter', 'latex');
legend([h1, h2], {'Mean', 'Mean $\pm$ Std'}, 'Interpreter', 'latex',   'Location', 'NW');
