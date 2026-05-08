%% SFVM and ROM for 1D Euler equations with two stochastic variables
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

% Ny_list = [4;8;16;32]
% state_time_list = [];
% state_step_list = [];
% flux_time_list = [];
% flux_step_list = [];
% difference_list = [];
% 
% for r = 1:length(Ny_list)
% Ny = Ny_list(r);

Ny = 16;
xpts = linspace(0,1,Nx+1);
ypts = linspace(0,1,Ny+1);
zpts = linspace(0,1,Ny+1);

x = 0.5*(xpts(1:Nx)+xpts(2:Nx+1));
y = 0.5*(ypts(1:Ny)+ypts(2:Ny+1));
z = 0.5*(zpts(1:Ny)+zpts(2:Ny+1));

[Y,X,Z] = meshgrid(y,x,z);


t0=0;
tmax=0.2;
tspan = [t0 tmax];

u = zeros(Nx,Ny,Ny);
rho = ones(Nx,Ny,Ny) * 0.125;
p = ones(Nx,Ny,Ny) * 0.1;


rho(X < 0.475 + 0.05*Y) = 1; 
p(X < 0.475 + 0.05*Y) = 1;  

E = zeros(Nx,Ny,Ny);
gamma = zeros(Ny,1);
for k = 1:Ny
    gamma(k) = 1.2 + 0.4 * y(k);
    E(:,:,k) = p(:,:,k) / (gamma(k) - 1) + 0.5 * rho(:,:,k) .* u(:,:,k).^2;
end

gamma_3D = reshape(gamma,1,1,Ny);

ics = [rho; rho .* u; E];

dx = (xpts(end)-xpts(1))/Nx;
dy = 1/Ny;


% Assemble parameters
params.Ny = Ny;
params.Nx = Nx;
params.dx = dx;
params.dy = dy;
params.gamma = gamma_3D;

% SFV - WENO with reconstructed states
options = odeset('RelTol',1e-6,'AbsTol',1e-8);

% tic
[~,U_state] = ode45(@(t,U) rhs_3D_euler_state(t,U,params), tspan, ics, options);
% toc

% f1 = @() ode45(@(t,U) rhs_3D_euler_state(t,U,params), tspan, ics, options);
% state_time_list = [state_time_list;timeit(f1)];
% state_step_list = [state_step_list;size(U_state,1)];

sol_state =reshape(U_state(end,:),3 * Nx,Ny,Ny);
sol_state_rho = sol_state(1:Nx,:,:);
sol_state_rhou = sol_state(Nx+1:2*Nx,:,:);
sol_state_u = sol_state_rhou ./ sol_state_rho;

tic
[~,U_flux] = ode45(@(t,U) rhs_3D_euler_flux(t,U,params), tspan, ics, options);
toc

% f2 = @() ode45(@(t,U) rhs_3D_euler_flux(t,U,params), tspan, ics, options);
% flux_time_list = [flux_time_list;timeit(f2)];
% flux_step_list = [flux_step_list;size(U_flux,1)];


sol_flux =reshape(U_flux(end,:),3 * Nx,Ny,Ny);
sol_flux_rho = sol_flux(1:Nx,:,:);
sol_flux_rhou = sol_flux(Nx+1:2*Nx,:,:);
sol_flux_u = sol_flux_rhou ./ sol_flux_rho;

er_flux = sum(abs(sol_state - sol_flux)...
    * (dx*dy^2),"all") / sum(abs(sol_state) * (dx*dy^2),"all");

% difference_list = [difference_list;er_flux];
disp("The relative L1 difference (recosntructed states and reconstructed fluxes) is " + er_flux)
% end

% Compute and plot mean and std of solution (with reconstructed fluxes)
sol_mean = sum(sum(sol_flux_rho, 2), 3) * dy^2;  
sol_diff = sol_flux_rho - sol_mean;  
sol_var = sum(sum(sol_diff.^2, 2), 3) / (Ny^2 - 1); 
sol_std = sqrt(sol_var); 
% 
figure
h1 = plot(x, sol_mean, 'LineWidth', 2, 'Color', 'r'); 
hold on
h2 = plot(x, sol_mean + sol_std, 'LineWidth', 2, 'Color', 'r', 'LineStyle', '--');
plot(x, sol_mean - sol_std, 'LineWidth', 2, 'Color', 'r', 'LineStyle', '--');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$\rho$', 'Interpreter', 'latex');
legend([h1, h2], {'Mean', 'Mean $\pm$ Std'}, 'Interpreter', 'latex',   'Location', 'best');


sol_mean = sum(sum(sol_flux_u, 2), 3) * dy^2;  
sol_diff = sol_flux_u - sol_mean;  
sol_var = sum(sum(sol_diff.^2, 2), 3) / (Ny^2 - 1); 
sol_std = sqrt(sol_var); 
% 
figure
h1 = plot(x, sol_mean, 'LineWidth', 2, 'Color', 'r'); 
hold on
h2 = plot(x, sol_mean + sol_std, 'LineWidth', 2, 'Color', 'r', 'LineStyle', '--');
plot(x, sol_mean - sol_std, 'LineWidth', 2, 'Color', 'r', 'LineStyle', '--');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$u$', 'Interpreter', 'latex');
legend([h1, h2], {'Mean', 'Mean $\pm$ Std'}, 'Interpreter', 'latex',   'Location', 'best');



%% ROM - POD
% 1D sampling
tspan_sample = linspace(t0, tmax, 100);
F_sample = zeros(Ny, Ny, 2*3*Nx*length(tspan_sample));
for i = 1:length(y)
    for j = 1:length(z)
        u = zeros(Nx,1);
        rho = ones(Nx,1) * 0.125;
        p = ones(Nx,1) * 0.1;
        rho(x < 0.475+ 0.05*y(i)) = 1; 
        p(x < 0.475+ 0.05*y(i)) = 1;     
        gamma_sampling = gamma(j);
        E = p / (gamma_sampling - 1) + 0.5 * rho .* u.^2;
        ics_sample = [rho; rho.*u; E];
        F_snap =  sampling_euler(Nx, gamma_sampling, ics_sample, tspan_sample);
        F_sample(i,j,:)= F_snap;
    end
end

a1 = 0.5*(-1/sqrt(3)); 
a2 = 0.5*(1/sqrt(3));

Fq_y1 = WENO_pages(F_sample,a1,"outflow");
Fq_y2 = WENO_pages(F_sample,a2,"outflow");

Fq_1 = pagetranspose(WENO_pages(pagetranspose(Fq_y1),a1,"outflow"));
Fq_2 = pagetranspose(WENO_pages(pagetranspose(Fq_y1),a2,"outflow"));
Fq_3 = pagetranspose(WENO_pages(pagetranspose(Fq_y2),a1,"outflow"));
Fq_4 = pagetranspose(WENO_pages(pagetranspose(Fq_y2),a2,"outflow"));
F_sample_q = zeros(4*Ny^2, size(F_sample,3));
F_sample_q(1:4:end, :) = reshape(Fq_1,Ny^2,[]); 
F_sample_q(2:4:end, :) = reshape(Fq_2,Ny^2,[]);
F_sample_q(3:4:end, :) = reshape(Fq_3,Ny^2,[]);
F_sample_q(4:4:end, :) = reshape(Fq_4,Ny^2,[]);

[Un,s,~] = svd(F_sample_q,"econ");

% ROM_er_list = [];
% ROM_time_list = [];
% ROM_step_list = [];
% N_list = 15:5:40;
% % % 
% for m = 1:length(N_list)
% 
% Nmode = N_list(m);
Nmode = 30;
Vn = Un(:,1:Nmode);

Bn = zeros(Ny^2,Nmode);
for i = 1:Ny^2
    Bn(i,:) = (Vn(4*i-3,:) +Vn(4*i-2,:)+Vn(4*i-1,:)+ Vn(4*i,:)) /4;
end

% Q-DEIM for hyper-reduction
[~,~,P] = qr(Vn',"vector");
NHR = Nmode;
ids = P(1:NHR);
Vn_HR = Vn(ids,:);


[y_ids, z_ids, q12_y_ids, q34_y_ids, q13_z_ids, q24_z_ids, map_q] = HR_ids_2s(ids,Ny);
params.gamma = gamma_3D(:,:,z_ids);

%
params.Vinv = pinv(Vn_HR);
params.B = Bn;
params.ids = ids;
params.y_ids = y_ids;
params.z_ids = z_ids;
params.q12_y_ids = q12_y_ids;
params.q34_y_ids = q34_y_ids;
params.q13_z_ids = q13_z_ids;
params.q24_z_ids = q24_z_ids;
params.map_q = map_q;

% tic
[~,U_ROM] = ode45(@(t,U) rhs_3D_euler_interp_HR(t,U,params), tspan, ics, options);
% toc

% 
% f3 = @() ode45(@(t,U) rhs_3D_euler_interp_HR(t,U,params), tspan, ics, options);
% ROM_time_list = [ROM_time_list; timeit(f3)];
% ROM_step_list = [ROM_step_list; size(U_ROM,1)];


sol_ROM = reshape(U_ROM(end,:),3*Nx,Ny,Ny);
sol_ROM_rho = sol_ROM(1:Nx,:,:);
sol_ROM_rhou = sol_ROM(Nx+1:2*Nx,:,:);
sol_ROM_u = sol_ROM_rhou ./ sol_ROM_rho;

% Compute ROM relative error
er_ROM = sum(abs(sol_flux - sol_ROM) * (dx*dy),"all") ...
    / sum(abs(sol_flux) * (dx*dy),"all");

% ROM_er_list = [ROM_er_list; er_ROM];
% end

disp("The relative L1 error (reconstructed fluxes and ROM) is " + er_ROM)


[YY,XX] = meshgrid(y,x);
figure
surf(XX,YY,reshape(sol_flux_rho(:,:,end),[],Ny))
xlabel('$x$'), ylabel('$y_1$'), zlabel('$\rho$')
zlim([-0.15 1.45])
view(52.5, 30)


% Singular value plot
figure
semilogy(diag(s), 'o-');
xlabel('Mode index');
ylabel('Singular value');
% xlim([0 100])

% ROM error and runtime plots
figure
semilogy(N_list, ROM_er_list, 'o-')
xlabel('Mode index')
ylabel('Error')

figure
scale = 1e3;
plot(N_list, scale * ROM_time_list./ROM_step_list, 'o-')
hold on
flux_time_Ny16 = flux_time_list(3)/flux_step_list(3);
plot(N_list, scale * flux_time_Ny16 * ones(length(N_list)),'r')
xlabel('Mode index')
ylabel('Runtime per step (ms)')
legend('ROM','Full SFV','Location','NW')


figure
plot(N_list, time_list, 'o-')
xlabel('Mode index');
ylabel('Runtime');
xlim([15 40])
ylim([0 250])



% Compute and plot mean and std of ROM solutions
sol_mean = sum(sum(sol_flux_rho, 2), 3) * dy^2;  
sol_diff = sol_flux_rho - sol_mean;  
sol_var = sum(sum(sol_diff.^2, 2), 3) / (Ny^2 - 1); 
sol_std = sqrt(sol_var); 
% 
figure
h1 = plot(x, sol_mean, 'LineWidth', 2, 'Color', 'r'); 
hold on
h2 = plot(x, sol_mean + sol_std, 'LineWidth', 2, 'Color', 'r', 'LineStyle', '--');
plot(x, sol_mean - sol_std, 'LineWidth', 2, 'Color', 'r', 'LineStyle', '--');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$\rho$', 'Interpreter', 'latex');
legend([h1, h2], {'Mean', 'Mean $\pm$ Std'}, 'Interpreter', 'latex',   'Location', 'best');


sol_mean = sum(sum(sol_flux_u, 2), 3) * dy^2;  
sol_diff = sol_flux_u - sol_mean;  
sol_var = sum(sum(sol_diff.^2, 2), 3) / (Ny^2 - 1); 
sol_std = sqrt(sol_var); 
% 
figure
h1 = plot(x, sol_mean, 'LineWidth', 2, 'Color', 'r'); 
hold on
h2 = plot(x, sol_mean + sol_std, 'LineWidth', 2, 'Color', 'r', 'LineStyle', '--');
plot(x, sol_mean - sol_std, 'LineWidth', 2, 'Color', 'r', 'LineStyle', '--');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$u$', 'Interpreter', 'latex');
legend([h1, h2], {'Mean', 'Mean $\pm$ Std'}, 'Interpreter', 'latex',   'Location', 'NW');


