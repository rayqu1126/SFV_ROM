%% SFVM and ROM for 1D Euler equations with one stochastic variables
clear
close all
addpath('func')

% Set default plot
set(0, 'defaultaxesfontsize',24,'defaultaxeslinewidth',2,...
       'defaultlinelinewidth',2,'defaultpatchlinewidth',2,...
       'defaulttextfontsize',24,'defaulttextinterpreter','latex');
set(0, 'DefaultFigurePosition',[0 0 600 400]);

% Initialize mesh (Nx: physical Ny: stochastic)
Nx = 256;

% Ny_list = [4;8;16;32]
% state_time_list = [];
% state_step_list = [];
% flux_time_list = [];
% flux_step_list = [];
% difference_list = [];
% 
% for r = 1:length(Ny_list)
% Ny = Ny_list(r);

Ny = 32;

% xpts = linspace(0,1,Nx+1);
xpts = linspace(-5,5,Nx+1);
ypts = linspace(0,1,Ny+1);

x = 0.5*(xpts(1:Nx)+xpts(2:Nx+1));
y = 0.5*(ypts(1:Ny)+ypts(2:Ny+1));

[Y,X] = meshgrid(y,x);

% Initialize solution array
t0=0;

tmax = 0.2; % for Sod 
% tmax= 1.8; % for Shu-Osher
tspan = [t0 tmax];

% IC for Sod shock tube 
u = zeros(Nx,Ny);
rho = ones(Nx,Ny) * 0.125;
p = ones(Nx,Ny) * 0.1;
rho(X < 0.475 + 0.05*Y) = 1; 
p(X < 0.475 + 0.05*Y) = 1;   
% 
% IC for Shu-Osher shock tube
% u = ones(Nx,Ny) * 2.629369;
% rho = ones(Nx,Ny) * 3.857143;
% p = ones(Nx,Ny) * 10.3333;
% 
% u(X >= -4.1 + 0.2*Y) = 0;
% rho(X >= -4.1 + 0.2 * Y) = 1 + 0.2 * sin(5 * X(X >= -4.1 + 0.2 * Y));
% p(X >= -4.1 + 0.2*Y) = 1;   


gamma = 1.4;
E = p / (gamma - 1) + 0.5 * rho .* u.^2;

ics = [rho; rho .* u; E];


% Assemble parameters
dx = (xpts(end)-xpts(1))/Nx;
dy = 1/Ny;

params.Ny = Ny;
params.Nx = Nx;
params.dx = dx;
params.dy = dy;
params.gamma = gamma;

% SFV - with reconstructed states
options = odeset('NonNegative', numel(ics), 'RelTol',1e-6,'AbsTol',1e-8);
% tic
[~,U_state] = ode45(@(t,U) rhs_2D_euler_state(t,U,params), tspan, ics, options);
% toc

sol_state = reshape(U_state(end,:),3*Nx,Ny);
sol_state_rho = sol_state(1:Nx,:);

% f1 = @() ode45(@(t,U) rhs_2D_euler_state(t,U,params), tspan, ics, options);
% state_time_list = [state_time_list;timeit(f1)];
% state_step_list = [state_step_list;size(U_state,1)];

% SFV - with reconstructed fluxes
% tic
[~,U_flux] = ode45(@(t,U) rhs_2D_euler_flux(t,U,params), tspan, ics, options);
% toc

% f2 = @() ode45(@(t,U) rhs_2D_euler_flux(t,U,params), tspan, ics, options);
% flux_time_list = [flux_time_list;timeit(f2)];
% flux_step_list = [flux_step_list;size(U_flux,1)];

sol_flux = reshape(U_flux(end,:),3*Nx,Ny);
sol_flux_rho = sol_flux(1:Nx,:);
sol_flux_rhou = sol_flux(Nx+1:2*Nx,:);
sol_flux_u = sol_flux_rhou ./ sol_flux_rho;

% Compute relative difference
er_flux = sum(abs(sol_state - sol_flux) * (dx*dy),"all") ...
    / sum(abs(sol_state) * (dx*dy),"all");

% difference_list = [difference_list;er_flux];
disp("The relative L1 difference (recosntructed states and reconstructed fluxes) is " + er_flux)
% 
% end

figure
surf(X,Y,sol_state_rho)
xlabel('$x$'), ylabel('$y$'), zlabel('$\rho$')
% zlim([-0.15 1.45])
view(52.5, 30)
% 
% Solution plots (with reconstructed fluxes)
figure
surf(X,Y,sol_flux_rho)
xlabel('$x$'), ylabel('$y$'), zlabel('$\rho$')
% zlim([-0.15 1.45])
view(52.5, 30)

figure
surf(X,Y,sol_flux_u)
xlabel('$x$'), ylabel('$y$'), zlabel('$u$')
% zlim([-0.15 1.45])
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

F_sample = zeros(Ny, 2*3*Nx*length(tspan_sample));
for i = 1:length(y)
    % Sampling IC for Sod shock tube
    u = zeros(Nx,1);
    rho = ones(Nx,1) * 0.125;
    p = ones(Nx,1) * 0.1;
    rho(x < 0.475+ 0.05*y(i)) = 1; 
    p(x < 0.475+ 0.05*y(i)) = 1;     
    % 
    % Sampling IC for Shu-Osher shock tube
    % u = ones(Nx,1) * 2.629369;
    % rho = ones(Nx,1) * 3.857143;
    % p = ones(Nx,1) * 10.3333;
    % u(x >= -4.1 + 0.2*y(i)) = 0;
    % rho(x >= -4.1 + 0.2 * y(i)) = 1 + 0.2 * sin(5 * x(x >= -4.1 + 0.2 * y(i)));
    % p(x >= -4.1 + 0.2*y(i)) = 1;   

    E = p / (gamma - 1) + 0.5 * rho .* u.^2;
    ics_sample = [rho; rho.*u; E];
    F_sample(i,:)= sampling_euler(Nx, gamma, ics_sample, tspan_sample);
end
  a1 = 0.5*(-1/sqrt(3)); 
  a2 = 0.5*(1/sqrt(3));
  Fq1 = WENO_2Darray(F_sample,a1,"outflow");
  Fq2 = WENO_2Darray(F_sample,a2,"outflow");
  F_sample_q = zeros(2*Ny, size(F_sample,2));

F_sample_q(1:2:end, :) = Fq1; 
F_sample_q(2:2:end, :) = Fq2; 

% SVD 
[Un,s,~] = svd(F_sample_q,"econ");
% 
% ROM_er_list = [];
% ROM_time_list = [];
% ROM_step_list = [];
% N_list = 3:13;

% for m = 1:length(N_list)
% Nmode = N_list(m);
Nmode = 6;
Vn = Un(:,1:Nmode);

% % Integral of basis
Bn = zeros(Ny,Nmode);
for i = 1:Ny
    Bn(i,:) = (Vn(2*i-1,:) + Vn(2*i,:)) /2;
end

% Q-DEIM hyper-reduction
[~,~,P] = qr(Vn',"vector");

NHR = Nmode;
ids = P(1:NHR);

% ids = 1:2*Ny;

Vn_HR = Vn(ids,:);

% HR indexing
[y_ids,q1_ids,q2_ids,type_merged,perm] = HR_ids_1s(ids,Ny);


params.Vinv = pinv(Vn_HR);
params.B = Bn;
params.ids = ids;
params.y_ids = y_ids;
params.q1_ids = q1_ids;
params.q2_ids = q2_ids;
params.type_merged = type_merged;
params.perm = perm;

% tic
[~,U_ROM] = ode45(@(t,U) rhs_2D_euler_interp_HR(t,U,params), tspan, ics, options);
% toc

% 
% f3 = @() ode45(@(t,U) rhs_2D_euler_interp_HR(t,U,params), tspan, ics, options);
% ROM_time_list = [ROM_time_list; timeit(f3)];
% ROM_step_list = [ROM_step_list; size(U_ROM,1)];

sol_ROM = reshape(U_ROM(end,:),3*Nx,Ny);
sol_ROM_rho = sol_ROM(1:Nx,:);
sol_ROM_rhou = sol_ROM(Nx+1:2*Nx,:);
sol_ROM_u = sol_ROM_rhou ./ sol_ROM_rho;

% Compute ROM relative error
er_ROM = sum(abs(sol_flux - sol_ROM) * (dx*dy),"all") ...
    / sum(abs(sol_flux) * (dx*dy),"all");
% 
% ROM_er_list = [ROM_er_list; er_ROM];
% end


disp("The relative L1 error (reconstructed fluxes and ROM) is " + er_ROM)


% Singular value plot
figure
semilogy(diag(s), 'o-')
xlabel('Mode index')
ylabel('Singular value')

% ROM error and runtime plots
% figure
% scale = 1e3;
% plot(N_list, scale * ROM_time_list./ROM_step_list, 'o-')
% hold on
% flux_time_Ny32 = flux_time_list(4)/flux_step_list(4);
% plot(N_list, scale * flux_time_Ny32 * ones(length(N_list)),'r')
% xlabel('Mode index')
% ylabel('Runtime per step (ms)')
% legend('ROM','Full SFV','Location','NW')
% xlim([3 13])
% 
% figure
% semilogy(N_list, ROM_er_list, 'o-')
% xlabel('Mode index')
% ylabel('Error')
% xlim([3 13])
% ylim([1e-12 1e-2])


% ROM solution plots
figure
surf(X,Y,sol_ROM_rho)
xlabel('$x$'), ylabel('$y$'), zlabel('$\rho$')
% zlim([-0.15 4.15])
view(52.5, 30)

figure
surf(X,Y,sol_ROM_u)
xlabel('$x$'), ylabel('$y$'), zlabel('$u$')
% zlim([-0.15 3])
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
legend([h1, h2], {'Mean', 'Mean $\pm$ Std'}, 'Interpreter', 'latex',   'Location', 'best');
% 

