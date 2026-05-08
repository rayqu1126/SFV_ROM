%% 1D Burgers' equations x 2D stochastic variable (interpolation method)
function dUdt = rhs_3D_euler_interp_HR(~,U,params)
  Ny = params.Ny;
  Nx = params.Nx;
  dx = params.dx;
  B = params.B;
  Vinv = params.Vinv;
  ids = params.ids;
  gamma = params.gamma;
  y_ids = params.y_ids;
  z_ids = params.z_ids;
  q12_y_ids = params.q12_y_ids;
  q34_y_ids = params.q34_y_ids;
  q13_z_ids = params.q13_z_ids;
  q24_z_ids = params.q24_z_ids;
  map_q = params.map_q;

  U = reshape(U,3*Nx,Ny,Ny);
    rho = U(1:Nx,y_ids,z_ids);
    rhou = U(Nx+1:2*Nx, y_ids,z_ids);
    E = U(2*Nx+1:end,y_ids,z_ids);

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

    F_rho_im12_1 = pagetranspose(WENO_pages_HR(pagetranspose(F_rho_im12),a1,"outflow",q12_y_ids));
    F_rho_im12_2 = pagetranspose(WENO_pages_HR(pagetranspose(F_rho_im12),a2,"outflow",q34_y_ids));
    F_rhou_im12_1 = pagetranspose(WENO_pages_HR(pagetranspose(F_rhou_im12),a1,"outflow",q12_y_ids));
    F_rhou_im12_2 = pagetranspose(WENO_pages_HR(pagetranspose(F_rhou_im12),a2,"outflow",q34_y_ids));
    F_E_im12_1 = pagetranspose(WENO_pages_HR(pagetranspose(F_E_im12),a1,"outflow",q12_y_ids));
    F_E_im12_2 = pagetranspose(WENO_pages_HR(pagetranspose(F_E_im12),a2,"outflow",q34_y_ids));


    F_rho_ip12_1 = pagetranspose(WENO_pages_HR(pagetranspose(F_rho_ip12),a1,"outflow",q12_y_ids));
    F_rho_ip12_2 = pagetranspose(WENO_pages_HR(pagetranspose(F_rho_ip12),a2,"outflow",q34_y_ids));
    F_rhou_ip12_1 = pagetranspose(WENO_pages_HR(pagetranspose(F_rhou_ip12),a1,"outflow",q12_y_ids));
    F_rhou_ip12_2 = pagetranspose(WENO_pages_HR(pagetranspose(F_rhou_ip12),a2,"outflow",q34_y_ids));
    F_E_ip12_1 = pagetranspose(WENO_pages_HR(pagetranspose(F_E_ip12),a1,"outflow",q12_y_ids));
    F_E_ip12_2 = pagetranspose(WENO_pages_HR(pagetranspose(F_E_ip12),a2,"outflow",q34_y_ids));
 


    F_rho_im12_1_perm = permute(F_rho_im12_1,[3 2 1]);
    F_rho_im12_q1 = permute(WENO_pages_HR(F_rho_im12_1_perm,a1,"outflow",q13_z_ids),[3 2 1]);
    F_rho_im12_q2 = permute(WENO_pages_HR(F_rho_im12_1_perm,a2,"outflow",q24_z_ids),[3 2 1]);
    F_rho_im12_2_perm = permute(F_rho_im12_2,[3 2 1]);
    F_rho_im12_q3 = permute(WENO_pages_HR(F_rho_im12_2_perm,a1,"outflow",q13_z_ids),[3 2 1]);
    F_rho_im12_q4 = permute(WENO_pages_HR(F_rho_im12_2_perm,a2,"outflow",q24_z_ids),[3 2 1]);
    
    F_rhou_im12_1_perm = permute(F_rhou_im12_1,[3 2 1]);
    F_rhou_im12_q1 = permute(WENO_pages_HR(F_rhou_im12_1_perm,a1,"outflow",q13_z_ids),[3 2 1]);
    F_rhou_im12_q2 = permute(WENO_pages_HR(F_rhou_im12_1_perm,a2,"outflow",q24_z_ids),[3 2 1]);
    F_rhou_im12_2_perm = permute(F_rhou_im12_2,[3 2 1]);
    F_rhou_im12_q3 = permute(WENO_pages_HR(F_rhou_im12_2_perm,a1,"outflow",q13_z_ids),[3 2 1]);
    F_rhou_im12_q4 = permute(WENO_pages_HR(F_rhou_im12_2_perm,a2,"outflow",q24_z_ids),[3 2 1]);
    

    F_E_im12_1_perm = permute(F_E_im12_1,[3 2 1]);
    F_E_im12_q1 = permute(WENO_pages_HR(F_E_im12_1_perm,a1,"outflow",q13_z_ids),[3 2 1]);
    F_E_im12_q2 = permute(WENO_pages_HR(F_E_im12_1_perm,a2,"outflow",q24_z_ids),[3 2 1]);
    F_E_im12_2_perm = permute(F_E_im12_2,[3 2 1]);
    F_E_im12_q3 = permute(WENO_pages_HR(F_E_im12_2_perm,a1,"outflow",q13_z_ids),[3 2 1]);
    F_E_im12_q4 = permute(WENO_pages_HR(F_E_im12_2_perm,a2,"outflow",q24_z_ids),[3 2 1]);
    
     
    Fim12_combined = zeros(length(ids), 3*Nx);
    for qtype = 1:4
        switch qtype
            case 1
                y_set = q12_y_ids;
                z_set = q13_z_ids;
                F_rho = F_rho_im12_q1;
                F_rhou = F_rhou_im12_q1;
                F_E = F_E_im12_q1;
    
            case 2
                y_set = q12_y_ids;
                z_set = q24_z_ids;
                F_rho = F_rho_im12_q2;
                F_rhou = F_rhou_im12_q2;
                F_E = F_E_im12_q2;
    
            case 3
                y_set = q34_y_ids;
                z_set = q13_z_ids;
                F_rho = F_rho_im12_q3;
                F_rhou = F_rhou_im12_q3;
                F_E = F_E_im12_q3;
    
            case 4
                y_set = q34_y_ids;
                z_set = q24_z_ids;
                F_rho = F_rho_im12_q4;
                F_rhou = F_rhou_im12_q4;
                F_E = F_E_im12_q4;
        end
    
        % loop over reduced tensor
        for iy = 1:length(y_set)
            y_global = y_ids(y_set(iy));
    
            for jz = 1:length(z_set)
                z_global = z_ids(z_set(jz));
    
                cell_lin = sub2ind([Ny,Ny], y_global, z_global);
                k = map_q(cell_lin, qtype);
    
                if k ~= 0
                    Fim12_combined(k,:) = [
                        F_rho(:,iy,jz);
                        F_rhou(:,iy,jz);
                        F_E(:,iy,jz)
                    ]';
                end
            end
        end
    end


    
    F_rho_ip12_1_perm = permute(F_rho_ip12_1,[3 2 1]);
    F_rho_ip12_q1 = permute(WENO_pages_HR(F_rho_ip12_1_perm,a1,"outflow",q13_z_ids),[3 2 1]);
    F_rho_ip12_q2 = permute(WENO_pages_HR(F_rho_ip12_1_perm,a2,"outflow",q24_z_ids),[3 2 1]);
    F_rho_ip12_2_perm = permute(F_rho_ip12_2,[3 2 1]);
    F_rho_ip12_q3 = permute(WENO_pages_HR(F_rho_ip12_2_perm,a1,"outflow",q13_z_ids),[3 2 1]);
    F_rho_ip12_q4 = permute(WENO_pages_HR(F_rho_ip12_2_perm,a2,"outflow",q24_z_ids),[3 2 1]);
    
    F_rhou_ip12_1_perm = permute(F_rhou_ip12_1,[3 2 1]);
    F_rhou_ip12_q1 = permute(WENO_pages_HR(F_rhou_ip12_1_perm,a1,"outflow",q13_z_ids),[3 2 1]);
    F_rhou_ip12_q2 = permute(WENO_pages_HR(F_rhou_ip12_1_perm,a2,"outflow",q24_z_ids),[3 2 1]);
    F_rhou_ip12_2_perm = permute(F_rhou_ip12_2,[3 2 1]);
    F_rhou_ip12_q3 = permute(WENO_pages_HR(F_rhou_ip12_2_perm,a1,"outflow",q13_z_ids),[3 2 1]);
    F_rhou_ip12_q4 = permute(WENO_pages_HR(F_rhou_ip12_2_perm,a2,"outflow",q24_z_ids),[3 2 1]);
    
    F_E_ip12_1_perm = permute(F_E_ip12_1,[3 2 1]);
    F_E_ip12_q1 = permute(WENO_pages_HR(F_E_ip12_1_perm,a1,"outflow",q13_z_ids),[3 2 1]);
    F_E_ip12_q2 = permute(WENO_pages_HR(F_E_ip12_1_perm,a2,"outflow",q24_z_ids),[3 2 1]);
    F_E_ip12_2_perm = permute(F_E_ip12_2,[3 2 1]);
    F_E_ip12_q3 = permute(WENO_pages_HR(F_E_ip12_2_perm,a1,"outflow",q13_z_ids),[3 2 1]);
    F_E_ip12_q4 = permute(WENO_pages_HR(F_E_ip12_2_perm,a2,"outflow",q24_z_ids),[3 2 1]);
    
    Fip12_combined = zeros(length(ids), 3*Nx);

for qtype = 1:4
    switch qtype
        case 1
            y_set = q12_y_ids;
            z_set = q13_z_ids;
            F_rho = F_rho_ip12_q1;
            F_rhou = F_rhou_ip12_q1;
            F_E = F_E_ip12_q1;

        case 2
            y_set = q12_y_ids;
            z_set = q24_z_ids;
            F_rho = F_rho_ip12_q2;
            F_rhou = F_rhou_ip12_q2;
            F_E = F_E_ip12_q2;

        case 3
            y_set = q34_y_ids;
            z_set = q13_z_ids;
            F_rho = F_rho_ip12_q3;
            F_rhou = F_rhou_ip12_q3;
            F_E = F_E_ip12_q3;

        case 4
            y_set = q34_y_ids;
            z_set = q24_z_ids;
            F_rho = F_rho_ip12_q4;
            F_rhou = F_rhou_ip12_q4;
            F_E = F_E_ip12_q4;
    end

    % loop over reduced tensor
    for iy = 1:length(y_set)
        y_global = y_ids(y_set(iy));

        for jz = 1:length(z_set)
            z_global = z_ids(z_set(jz));

            cell_lin = sub2ind([Ny,Ny], y_global, z_global);
            k = map_q(cell_lin, qtype);

            if k ~= 0
                Fip12_combined(k,:) = [
                    F_rho(:,iy,jz);
                    F_rhou(:,iy,jz);
                    F_E(:,iy,jz)
                ]';
            end
        end
    end
end


   intFim12 = B * Vinv * Fim12_combined; intFim12 = intFim12';
   intFip12 = B * Vinv * Fip12_combined; intFip12 = intFip12';

   dUdt = reshape(-1/dx*(intFip12 - intFim12), [],1);

end
