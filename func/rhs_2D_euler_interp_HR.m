%% 1D Euler equations x 1D stochastic variable (interpolation method)
function dUdt = rhs_2D_euler_interp_HR(~,U,params)
  Ny = params.Ny;
  Nx = params.Nx;
  dx = params.dx;
  B = params.B;
  Vinv = params.Vinv;
  ids = params.ids;
  gamma = params.gamma;
  y_ids = params.y_ids;
  q1_ids = params.q1_ids;
  q2_ids = params.q2_ids;
  type_merged = params.type_merged;
  perm = params.perm;

 
  % Unpack variables
  U = reshape(U,3 * Nx,Ny);
  rho = U(1:Nx,y_ids);
  rhou = U(Nx+1:2*Nx, y_ids);
  E = U(2*Nx+1:end,y_ids);

  u = rhou ./ rho;
  p = (gamma - 1) * (E - 0.5 * rho .* u.^2);



  
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
  F_rho_im12_1 = WENO_2Darray_HR(F_rho_im12',a1,"outflow",q1_ids)';
  F_rho_im12_2 = WENO_2Darray_HR(F_rho_im12',a2,"outflow",q2_ids)';
  F_rhou_im12_1 = WENO_2Darray_HR(F_rhou_im12',a1,"outflow",q1_ids)';
  F_rhou_im12_2 = WENO_2Darray_HR(F_rhou_im12',a2,"outflow",q2_ids)';
  F_E_im12_1 = WENO_2Darray_HR(F_E_im12',a1,"outflow",q1_ids)';
  F_E_im12_2 = WENO_2Darray_HR(F_E_im12',a2,"outflow",q2_ids)';


  Fim12_q1 = [F_rho_im12_1; F_rhou_im12_1; F_E_im12_1];
  Fim12_q2 = [F_rho_im12_2; F_rhou_im12_2; F_E_im12_2];

  Fim12_combined = zeros(length(ids),3*Nx);
  i1 = 1; i2 = 1;
  for k = 1:length(type_merged)
    if type_merged(k) == 1
        Fim12_combined(k,:) = Fim12_q1(:,i1);
        i1 = i1 + 1;
    else
        Fim12_combined(k,:) = Fim12_q2(:,i2);
        i2 = i2 + 1;
    end
  end

  F_rho_ip12_1 = WENO_2Darray_HR(F_rho_ip12',a1,"outflow",q1_ids)';
  F_rho_ip12_2 = WENO_2Darray_HR(F_rho_ip12',a2,"outflow",q2_ids)';
  F_rhou_ip12_1 = WENO_2Darray_HR(F_rhou_ip12',a1,"outflow",q1_ids)';
  F_rhou_ip12_2 = WENO_2Darray_HR(F_rhou_ip12',a2,"outflow",q2_ids)';
  F_E_ip12_1 = WENO_2Darray_HR(F_E_ip12',a1,"outflow",q1_ids)';
  F_E_ip12_2 = WENO_2Darray_HR(F_E_ip12',a2,"outflow",q2_ids)';

  Fip12_q1 = [F_rho_ip12_1; F_rhou_ip12_1; F_E_ip12_1];
  Fip12_q2 = [F_rho_ip12_2; F_rhou_ip12_2; F_E_ip12_2];
    
  % Combine flux in order of quadrature
  Fip12_combined = zeros(length(ids),3*Nx);
  i1 = 1; i2 = 1;
  for k = 1:length(type_merged)
    if type_merged(k) == 1
        Fip12_combined(k,:) = Fip12_q1(:,i1);
        i1 = i1 + 1;
    else
        Fip12_combined(k,:) = Fip12_q2(:,i2);
        i2 = i2 + 1;
    end
  end



  % Flux integrals with reduced basis and HR indices
  intFim12 = B * Vinv * Fim12_combined(perm,:); intFim12 = intFim12';
  intFip12 = B * Vinv * Fip12_combined(perm,:); intFip12 = intFip12';

  dUdt = reshape(-1/dx*(intFip12 - intFim12), [],1);
end