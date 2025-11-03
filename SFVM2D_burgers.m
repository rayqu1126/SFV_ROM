%% SFVM and ROM for 1D Burgers equations with 1D stochstic variables
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
Ny = 16;

xpts = linspace(0,1,Nx+1);
ypts = linspace(0,1,Ny+1);

x = 0.5*(xpts(1:Nx)+xpts(2:Nx+1));
y = 0.5*(ypts(1:Ny)+ypts(2:Ny+1));

[Y,X] = meshgrid(y,x);

% Initialize solution array
ics = sin(2*pi*X)+ 0.5 * Y;

t0=0;
tmax=0.3;
tspan = [t0 tmax];

% Assemble parameters
dx = 1/Nx;
dy = 1/Ny;
params.N = Ny;
params.Nx = Nx;
params.dx = 1/Nx;
params.dy = 1/Ny;

% SFV - WENO with reconstructed states
options = odeset('RelTol',1e-6,'AbsTol',1e-8);

tic
[~,U_state] = ode45(@(t,U) rhs_2D_burgers_state(t,U,params), tspan, ics, options);
toc

sol_state = reshape(U_state(end,:),Nx,Ny);

figure
surf(X,Y,sol_state)
xlabel('X'), ylabel('Y'), zlabel('SFV (with reconstructed states)')

% SFVM - flux reconstruction
tic
[~,U_flux] = ode45(@(t,U) rhs_2D_burgers_flux(t,U,params), tspan, ics, options);
toc

sol_flux = reshape(U_flux(end,:),Nx,Ny);

figure
surf(X,Y,sol_flux)
xlabel('X'), ylabel('Y'), zlabel('SFV (with reconstructed fluxes)')

er_flux = sum(abs(sol_state - sol_flux)...
    * (dx*dy),"all") / sum(abs(sol_state) * (dx*dy),"all");

disp("The relative L1 difference (recosntructed states and fluxes) is " + er_flux)
%% ROM - POD
% 1D sampling
tspan_sample = linspace(t0, tmax, 200);
F_sample = zeros(2*Ny, Nx*length(tspan_sample));

quadL = y + 0.5*(-1/sqrt(3)) * dy;
quadR = y + 0.5*(1/sqrt(3)) * dy;   
quad_combined = reshape([quadL; quadR], [], 1)';
for i = 1:length(quad_combined)
    ics_sample = sin(2*pi*x) + 0.5 * quad_combined(i);
    F_sample(i,:)= sampling_burgers(Nx, ics_sample, tspan_sample);
end

% Construct reduced basis
Nmode = 30;
[Vn,s,~] = svd(F_sample);
Vn = Vn(:,1:Nmode);

% Integral of basis
Bn = zeros(Ny,Nmode);
for i = 1:Ny
    Bn(i,:) = (Vn(2*i-1,:) + Vn(2*i,:)) /2;
end

% Q-DEIM hyper-reduction
[~,~,P] = qr(Vn',"vector");

NHR = Nmode;
ids = P(1:NHR);
Vn_HR = Vn(ids,:);

%
params.Vinv = pinv(Vn_HR);
params.B = Bn;
params.ids = ids;

[~,U_ROM] = ode45(@(t,U) rhs_2D_burgers_interp(t,U,params), tspan, ics, options);
%

sol_ROM = reshape(U_ROM(end,:),Nx,Ny);

figure
surf(X,Y,sol_ROM)
xlabel('X'), ylabel('Y'), zlabel('SFV (ROM)')

er_ROM = sum(abs(sol_flux - sol_ROM)...
    * (dx*dy),"all") / sum(abs(sol_flux) * (dx*dy),"all");

disp("The relative L1 error (reconstructed fluxes and ROM) is " + er_ROM)
