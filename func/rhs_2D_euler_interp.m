%% 1D Euler equations x 1D stochastic variable (interpolation method)
function dUdt = rhs_2D_euler_interp(~,U,params)
  Ny = params.Ny;
  Nx = params.Nx;
  dx = params.dx;
  B = params.B;
  Vinv = params.Vinv;
  ids = params.ids;
  gamma = params.gamma;
  % Vn = params.Vn;

  % Unpack variables
  U = reshape(U,3 * Nx,Ny);
  rho = U(1:Nx,:);
  rhou = U(Nx+1:2*Nx, :);
  E = U(2*Nx+1:end,:);
    
  u = rhou ./ rho;
  p = (gamma - 1) * (E - 0.5 * rho .* u.^2);

  % Reconstruct states in physical domain
  % [rhoL,rhoR] = WENO_2D_reflect(rho,Nx,Ny,dx,"rho");
  % [uL,uR] = WENO_2D_reflect(u,Nx,Ny,dx,"u");
  % [pL,pR] = WENO_2D_reflect(p,Nx,Ny,dx,"p");



  rhoL = WENO_2Darray(rho,-0.5,"outflow"); 
  rhoR = WENO_2Darray(rho,0.5,"outflow");
  uL = WENO_2Darray(u,-0.5,"outflow"); 
  uR = WENO_2Darray(u,0.5,"outflow");
  pL = WENO_2Darray(p,-0.5,"outflow"); 
  pR = WENO_2Darray(p,0.5,"outflow");


  % Flux on left interface (reflective wall)
  rhoim1 = [rhoR(1,:);rhoR(1:Nx-1,:)];
  rhoi = rhoL;
  uim1 = [uR(1,:);uR(1:Nx-1,:)];
  ui = uL;
  pim1 = [pR(1,:);pR(1:Nx-1,:)];
  pi = pL;

  Eim1 = pim1 / (gamma - 1) + 0.5 * rhoim1 .* uim1.^2;
  Ei = pi / (gamma - 1) + 0.5 * rhoi .* ui.^2;    

  lambda = max(abs(ui)+real(sqrt(gamma*pi./rhoi)),...
        abs(uim1)+real(sqrt(gamma*pim1./rhoim1)));
    
  F_rho_im12 = (rhoi.*ui + rhoim1.*uim1)/2 - 0.5*lambda.*(rhoi-rhoim1);
  F_rhou_im12 = (rhoi.*ui.^2+pi + rhoim1.*uim1.^2+pim1)/2 - 0.5*lambda.*(rhoi .* ui- rhoim1 .* uim1);
  F_E_im12 = (ui.*(Ei + pi) + uim1.*(Eim1 + pim1))/2 - 0.5*lambda.*(Ei-Eim1);

  % Flux on right interface (reflective wall)
  rhoip1 = [rhoL(2:Nx,:);rhoL(Nx,:)];
  rhoi = rhoR;
  uip1 = [uL(2:Nx,:);uL(Nx,:)];
  ui = uR;
  pip1 = [pL(2:Nx,:);pL(Nx,:)];
  pi = pR;

  Eip1 = pip1 / (gamma - 1) + 0.5 * rhoip1 .* uip1.^2;
  Ei = pi / (gamma - 1) + 0.5 * rhoi .* ui.^2;    

   lambda = max(abs(ui)+real(sqrt(gamma*pi./rhoi)),...
        abs(uip1)+real(sqrt(gamma*pip1./rhoip1)));


  F_rho_ip12 = (rhoi.*ui + rhoip1.*uip1)/2 - 0.5*lambda.*(rhoip1-rhoi);
  F_rhou_ip12 = (rhoi.*ui.^2+pi + rhoip1.*uip1.^2+pip1)/2 - 0.5*lambda.*(rhoip1 .* uip1- rhoi .* ui);
  F_E_ip12 = (ui.*(Ei + pi) + uip1.*(Eip1 + pip1))/2 - 0.5*lambda.*(Eip1-Ei);
   
  % Reconstruct flux in quadrature points
  a1 = 0.5*(-1/sqrt(3)); 
  a2 = 0.5*(1/sqrt(3));
  % 
  % F_rho_im12_1 = (WENO_QP(Ny,a1)*F_rho_im12')';
  % F_rho_im12_2 = (WENO_QP(Ny,a2)*F_rho_im12')';
  % F_rhou_im12_1 = (WENO_QP(Ny,a1)*F_rhou_im12')';
  % F_rhou_im12_2 = (WENO_QP(Ny,a2)*F_rhou_im12')';
  % F_E_im12_1 = (WENO_QP(Ny,a1)*F_E_im12')';
  % F_E_im12_2 = (WENO_QP(Ny,a2)*F_E_im12')';


  F_rho_im12_1 = WENO_2Darray(F_rho_im12',a1,"outflow")';
  F_rho_im12_2 = WENO_2Darray(F_rho_im12',a2,"outflow")';
  F_rhou_im12_1 = WENO_2Darray(F_rhou_im12',a1,"outflow")';
  F_rhou_im12_2 = WENO_2Darray(F_rhou_im12',a2,"outflow")';
  F_E_im12_1 = WENO_2Darray(F_E_im12',a1,"outflow")';
  F_E_im12_2 = WENO_2Darray(F_E_im12',a2,"outflow")';


  Fim12_1 = [F_rho_im12_1; F_rhou_im12_1; F_E_im12_1];
  Fim12_2 = [F_rho_im12_2; F_rhou_im12_2; F_E_im12_2];

  Fim12_combined = zeros(2*Ny,3*Nx);
  for i = 1:3*Nx
      Fim12_combined(:,i) =  reshape([Fim12_1(i,:); Fim12_2(i,:)], [], 1)';
  end

  % 
  % F_rho_ip12_1 = (WENO_QP(Ny,a1)*F_rho_ip12')';
  % F_rho_ip12_2 = (WENO_QP(Ny,a2)*F_rho_ip12')';
  % F_rhou_ip12_1 = (WENO_QP(Ny,a1)*F_rhou_ip12')';
  % F_rhou_ip12_2 = (WENO_QP(Ny,a2)*F_rhou_ip12')';
  % F_E_ip12_1 = (WENO_QP(Ny,a1)*F_E_ip12')';
  % F_E_ip12_2 = (WENO_QP(Ny,a2)*F_E_ip12')';


  F_rho_ip12_1 = WENO_2Darray(F_rho_ip12',a1,"outflow")';
  F_rho_ip12_2 = WENO_2Darray(F_rho_ip12',a2,"outflow")';
  F_rhou_ip12_1 = WENO_2Darray(F_rhou_ip12',a1,"outflow")';
  F_rhou_ip12_2 = WENO_2Darray(F_rhou_ip12',a2,"outflow")';
  F_E_ip12_1 = WENO_2Darray(F_E_ip12',a1,"outflow")';
  F_E_ip12_2 = WENO_2Darray(F_E_ip12',a2,"outflow")';

  Fip12_1 = [F_rho_ip12_1; F_rhou_ip12_1; F_E_ip12_1];
  Fip12_2 = [F_rho_ip12_2; F_rhou_ip12_2; F_E_ip12_2];
    
  % Combine flux in order of quadrature
  Fip12_combined = zeros(2*Ny,3*Nx);
  for i = 1:3*Nx
    Fip12_combined(:,i) =  reshape([Fip12_1(i,:); Fip12_2(i,:)], [], 1)';
  end

  % Flux integrals with reduced basis and HR indices
  intFim12 = B * Vinv * Fim12_combined(ids,:); intFim12 = intFim12';
  intFip12 = B * Vinv * Fip12_combined(ids,:); intFip12 = intFip12';

  dUdt = reshape(-1/dx*(intFip12 - intFim12), [],1);
  
  % dUdt_flux = reshape(-1/dx*((Fip12_1+Fip12_2)/2-(Fim12_1+Fim12_2)/2), [],1);


  % 
  % rel_err_dF_ip12 = norm(-intFip12 + (Fip12_1+Fip12_2)/2) / norm((Fip12_1+Fip12_2)/2)
  % rel_err_dF_im12 = norm(-intFim12 + (Fim12_1+Fim12_2)/2) / norm((Fim12_1+Fim12_2)/2)
  % rel_err_dF = norm(dUdt - dUdt_flux) / norm(dUdt_flux)
  % rel_Fim12_Er = norm(Fim12_combined - Vn * Vn' * Fim12_combined) /  norm(Fim12_combined)
  % rel_Fip12_Er = norm(Fip12_combined - Vn * Vn' * Fip12_combined) / norm(Fip12_combined)
  % 
  % norm(norm(intFim12))
  % norm(norm(intFip12))
  % norm(norm(intFim12-intFip12))
end