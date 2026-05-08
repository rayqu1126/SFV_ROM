%% 1D Burgers' equations x 2D stochastic variable (interpolation method)
function dUdt = rhs_3D_euler_interp(~,U,params)
  Ny = params.Ny;
  Nx = params.Nx;
  dx = params.dx;
  B = params.B;
  Vinv = params.Vinv;
  ids = params.ids;
  gamma = params.gamma;


  U = reshape(U,3*Nx,Ny,Ny);
    rho = U(1:Nx,:,:);
    rhou = U(Nx+1:2*Nx, :,:);
    E = U(2*Nx+1:end,:,:);

    u = rhou ./ rho;
    p = (gamma - 1) .* (E - 0.5 * rho .* u.^2);

    rhoL = WENO_pages(rho,-0.5,"outflow"); 
    rhoR = WENO_pages(rho,0.5,"outflow");
    uL = WENO_pages(u,-0.5,"outflow"); 
    uR = WENO_pages(u,0.5,"outflow");
    pL = WENO_pages(p,-0.5,"outflow"); 
    pR = WENO_pages(p,0.5,"outflow");

    % Flux on left interface
    rhoim1 = [rhoR(1,:,:);rhoR(1:Nx-1,:,:)];
    rhoi = rhoL;
    uim1 = [uR(1,:,:);uR(1:Nx-1,:,:)];
    ui = uL;
    pim1 = [pR(1,:,:);pR(1:Nx-1,:,:)];
    pi = pL;

    Eim1 = pim1 ./ (gamma - 1) + 0.5 * rhoim1 .* uim1.^2;
    Ei = pi ./ (gamma - 1) + 0.5 * rhoi .* ui.^2;    

    lambda = max(abs(ui)+sqrt(gamma.*pi./rhoi),...
        abs(uim1)+sqrt(gamma.*pim1./rhoim1));
    
    F_rho_im12 = (rhoi.*ui + rhoim1.*uim1)/2 - 0.5*lambda.*(rhoi-rhoim1);
    F_rhou_im12 = (rhoi.*ui.^2+pi + rhoim1.*uim1.^2+pim1)/2 - 0.5*lambda.*(rhoi .* ui- rhoim1 .* uim1);
    F_E_im12 = (ui.*(Ei + pi) + uim1.*(Eim1 + pim1))/2 - 0.5*lambda.*(Ei-Eim1);

    % Flux on right interface
    rhoip1 = [rhoL(2:Nx,:,:);rhoL(Nx,:,:)];
    rhoi = rhoR;
    uip1 = [uL(2:Nx,:,:);uL(Nx,:,:)];
    ui = uR;
    pip1 = [pL(2:Nx,:,:);pL(Nx,:,:)];
    pi = pR;

    Eip1 = pip1 ./ (gamma - 1) + 0.5 * rhoip1 .* uip1.^2;
    Ei = pi ./ (gamma - 1) + 0.5 * rhoi .* ui.^2;    


    lambda = max(abs(ui)+sqrt(gamma.*pi./rhoi),...
        abs(uip1)+sqrt(gamma.*pip1./rhoip1));


    F_rho_ip12 = (rhoi.*ui + rhoip1.*uip1)/2 - 0.5*lambda.*(rhoip1-rhoi);
    F_rhou_ip12 = (rhoi.*ui.^2+pi + rhoip1.*uip1.^2+pip1)/2 - 0.5*lambda.*(rhoip1 .* uip1- rhoi .* ui);
    F_E_ip12 = (ui.*(Ei + pi) + uip1.*(Eip1 + pip1))/2 - 0.5*lambda.*(Eip1-Ei);
    
    
    % Reconstruct flux in quadrature points
    a1 = 0.5*(-1/sqrt(3)); 
    a2 = 0.5*(1/sqrt(3));
    F_rho_im12_1 = pagetranspose(WENO_pages(pagetranspose(F_rho_im12),a1,"outflow"));
    F_rho_im12_2 = pagetranspose(WENO_pages(pagetranspose(F_rho_im12),a2,"outflow"));
    F_rhou_im12_1 = pagetranspose(WENO_pages(pagetranspose(F_rhou_im12),a1,"outflow"));
    F_rhou_im12_2 = pagetranspose(WENO_pages(pagetranspose(F_rhou_im12),a2,"outflow"));
    F_E_im12_1 = pagetranspose(WENO_pages(pagetranspose(F_E_im12),a1,"outflow"));
    F_E_im12_2 = pagetranspose(WENO_pages(pagetranspose(F_E_im12),a2,"outflow"));


    F_rho_ip12_1 = pagetranspose(WENO_pages(pagetranspose(F_rho_ip12),a1,"outflow"));
    F_rho_ip12_2 = pagetranspose(WENO_pages(pagetranspose(F_rho_ip12),a2,"outflow"));
    F_rhou_ip12_1 = pagetranspose(WENO_pages(pagetranspose(F_rhou_ip12),a1,"outflow"));
    F_rhou_ip12_2 = pagetranspose(WENO_pages(pagetranspose(F_rhou_ip12),a2,"outflow"));
    F_E_ip12_1 = pagetranspose(WENO_pages(pagetranspose(F_E_ip12),a1,"outflow"));
    F_E_ip12_2 = pagetranspose(WENO_pages(pagetranspose(F_E_ip12),a2,"outflow"));


    F_rho_im12_1_perm = permute(F_rho_im12_1,[3 2 1]);
    F_rho_im12_q1 = permute(WENO_pages(F_rho_im12_1_perm,a1,"outflow"),[3 2 1]);
    F_rho_im12_q2 = permute(WENO_pages(F_rho_im12_1_perm,a2,"outflow"),[3 2 1]);
    F_rho_im12_2_perm = permute(F_rho_im12_2,[3 2 1]);
    F_rho_im12_q3 = permute(WENO_pages(F_rho_im12_2_perm,a1,"outflow"),[3 2 1]);
    F_rho_im12_q4 = permute(WENO_pages(F_rho_im12_2_perm,a2,"outflow"),[3 2 1]);

    F_rhou_im12_1_perm = permute(F_rhou_im12_1,[3 2 1]);
    F_rhou_im12_q1 = permute(WENO_pages(F_rhou_im12_1_perm,a1,"outflow"),[3 2 1]);
    F_rhou_im12_q2 = permute(WENO_pages(F_rhou_im12_1_perm,a2,"outflow"),[3 2 1]);
    F_rhou_im12_2_perm = permute(F_rhou_im12_2,[3 2 1]);
    F_rhou_im12_q3 = permute(WENO_pages(F_rhou_im12_2_perm,a1,"outflow"),[3 2 1]);
    F_rhou_im12_q4 = permute(WENO_pages(F_rhou_im12_2_perm,a2,"outflow"),[3 2 1]);

    F_E_im12_1_perm = permute(F_E_im12_1,[3 2 1]);
    F_E_im12_q1 = permute(WENO_pages(F_E_im12_1_perm,a1,"outflow"),[3 2 1]);
    F_E_im12_q2 = permute(WENO_pages(F_E_im12_1_perm,a2,"outflow"),[3 2 1]);
    F_E_im12_2_perm = permute(F_E_im12_2,[3 2 1]);
    F_E_im12_q3 = permute(WENO_pages(F_E_im12_2_perm,a1,"outflow"),[3 2 1]);
    F_E_im12_q4 = permute(WENO_pages(F_E_im12_2_perm,a2,"outflow"),[3 2 1]);


    F_rho_ip12_1_perm = permute(F_rho_ip12_1,[3 2 1]);
    F_rho_ip12_q1 = permute(WENO_pages(F_rho_ip12_1_perm,a1,"outflow"),[3 2 1]);
    F_rho_ip12_q2 = permute(WENO_pages(F_rho_ip12_1_perm,a2,"outflow"),[3 2 1]);
    F_rho_ip12_2_perm = permute(F_rho_ip12_2,[3 2 1]);
    F_rho_ip12_q3 = permute(WENO_pages(F_rho_ip12_2_perm,a1,"outflow"),[3 2 1]);
    F_rho_ip12_q4 = permute(WENO_pages(F_rho_ip12_2_perm,a2,"outflow"),[3 2 1]);

    F_rhou_ip12_1_perm = permute(F_rhou_ip12_1,[3 2 1]);
    F_rhou_ip12_q1 = permute(WENO_pages(F_rhou_ip12_1_perm,a1,"outflow"),[3 2 1]);
    F_rhou_ip12_q2 = permute(WENO_pages(F_rhou_ip12_1_perm,a2,"outflow"),[3 2 1]);
    F_rhou_ip12_2_perm = permute(F_rhou_ip12_2,[3 2 1]);
    F_rhou_ip12_q3 = permute(WENO_pages(F_rhou_ip12_2_perm,a1,"outflow"),[3 2 1]);
    F_rhou_ip12_q4 = permute(WENO_pages(F_rhou_ip12_2_perm,a2,"outflow"),[3 2 1]);

    F_E_ip12_1_perm = permute(F_E_ip12_1,[3 2 1]);
    F_E_ip12_q1 = permute(WENO_pages(F_E_ip12_1_perm,a1,"outflow"),[3 2 1]);
    F_E_ip12_q2 = permute(WENO_pages(F_E_ip12_1_perm,a2,"outflow"),[3 2 1]);
    F_E_ip12_2_perm = permute(F_E_ip12_2,[3 2 1]);
    F_E_ip12_q3 = permute(WENO_pages(F_E_ip12_2_perm,a1,"outflow"),[3 2 1]);
    F_E_ip12_q4 = permute(WENO_pages(F_E_ip12_2_perm,a2,"outflow"),[3 2 1]);


    Fim12_q1 = [F_rho_im12_q1; F_rhou_im12_q1; F_E_im12_q1];
    Fim12_q2 = [F_rho_im12_q2; F_rhou_im12_q2; F_E_im12_q2];
    Fim12_q3 = [F_rho_im12_q3; F_rhou_im12_q3; F_E_im12_q3];
    Fim12_q4 = [F_rho_im12_q4; F_rhou_im12_q4; F_E_im12_q4];


    Fip12_q1 = [F_rho_ip12_q1; F_rhou_ip12_q1; F_E_ip12_q1];
    Fip12_q2 = [F_rho_ip12_q2; F_rhou_ip12_q2; F_E_ip12_q2];
    Fip12_q3 = [F_rho_ip12_q3; F_rhou_ip12_q3; F_E_ip12_q3];
    Fip12_q4 = [F_rho_ip12_q4; F_rhou_ip12_q4; F_E_ip12_q4];


  Fip12_combined = zeros(4*Ny^2,3*Nx);
  Fip12_combined(1:4:end) = reshape(Fip12_q1,[],Ny^2)';
  Fip12_combined(2:4:end) = reshape(Fip12_q2,[],Ny^2)';
  Fip12_combined(3:4:end) = reshape(Fip12_q3,[],Ny^2)';
  Fip12_combined(4:4:end) = reshape(Fip12_q4,[],Ny^2)';


  Fim12_combined = zeros(4*Ny^2,3*Nx);
  Fim12_combined(1:4:end) = reshape(Fim12_q1,[],Ny^2)';
  Fim12_combined(2:4:end) = reshape(Fim12_q2,[],Ny^2)';
  Fim12_combined(3:4:end) = reshape(Fim12_q3,[],Ny^2)';
  Fim12_combined(4:4:end) = reshape(Fim12_q4,[],Ny^2)';

  intFim12 = B * Vinv * Fim12_combined(ids,:); intFim12 = intFim12';
  intFip12 = B * Vinv * Fip12_combined(ids,:); intFip12 = intFip12';

  dUdt = reshape(-1/dx*(intFip12 - intFim12), [],1);

end