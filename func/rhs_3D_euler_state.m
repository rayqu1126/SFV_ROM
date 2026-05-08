%% 1D Euler equations x 2D stochastic variable (state reconstruction)
function dUdt = rhs_3D_euler_state(~,U,params)
    Nx = params.Nx;
    Ny = params.Ny;
    dx = params.dx;
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


    a1 = 0.5*(-1/sqrt(3));
    a2 = 0.5*(1/sqrt(3));
    


        rhoL1 = pagetranspose(WENO_pages(pagetranspose(rhoL),a1,"outflow"));
        rhoL2 = pagetranspose(WENO_pages(pagetranspose(rhoL),a2,"outflow"));
        rhoR1 = pagetranspose(WENO_pages(pagetranspose(rhoR),a1,"outflow"));
        rhoR2 = pagetranspose(WENO_pages(pagetranspose(rhoR),a2,"outflow"));
        
        uL1 = pagetranspose(WENO_pages(pagetranspose(uL),a1,"outflow"));
        uL2 = pagetranspose(WENO_pages(pagetranspose(uL),a2,"outflow"));
        uR1 = pagetranspose(WENO_pages(pagetranspose(uR),a1,"outflow"));
        uR2 = pagetranspose(WENO_pages(pagetranspose(uR),a2,"outflow"));
        
        pL1 = pagetranspose(WENO_pages(pagetranspose(pL),a1,"outflow"));
        pL2 = pagetranspose(WENO_pages(pagetranspose(pL),a2,"outflow"));
        pR1 = pagetranspose(WENO_pages(pagetranspose(pR),a1,"outflow"));
        pR2 = pagetranspose(WENO_pages(pagetranspose(pR),a2,"outflow"));
        

        rhoL1_perm = permute(rhoL1,[3 2 1]); 
        rhoL11 = permute(WENO_pages(rhoL1_perm,a1,"outflow"),[3 2 1]);
        rhoL12 = permute(WENO_pages(rhoL1_perm,a2,"outflow"),[3 2 1]); 

        rhoL2_perm = permute(rhoL2,[3 2 1]); 
        rhoL21 = permute(WENO_pages(rhoL2_perm,a1,"outflow"),[3 2 1]);
        rhoL22 = permute(WENO_pages(rhoL2_perm,a2,"outflow"),[3 2 1]); 

        rhoR1_perm = permute(rhoR1,[3 2 1]); 
        rhoR11 = permute(WENO_pages(rhoR1_perm,a1,"outflow"),[3 2 1]);
        rhoR12 = permute(WENO_pages(rhoR1_perm,a2,"outflow"),[3 2 1]); 

        rhoR2_perm = permute(rhoR2,[3 2 1]); 
        rhoR21 = permute(WENO_pages(rhoR2_perm,a1,"outflow"),[3 2 1]);
        rhoR22 = permute(WENO_pages(rhoR2_perm,a2,"outflow"),[3 2 1]); 

        uL1_perm = permute(uL1,[3 2 1]); 
        uL11 = permute(WENO_pages(uL1_perm,a1,"outflow"),[3 2 1]);
        uL12 = permute(WENO_pages(uL1_perm,a2,"outflow"),[3 2 1]); 

        uL2_perm = permute(uL2,[3 2 1]); 
        uL21 = permute(WENO_pages(uL2_perm,a1,"outflow"),[3 2 1]);
        uL22 = permute(WENO_pages(uL2_perm,a2,"outflow"),[3 2 1]); 

        uR1_perm = permute(uR1,[3 2 1]); 
        uR11 = permute(WENO_pages(uR1_perm,a1,"outflow"),[3 2 1]);
        uR12 = permute(WENO_pages(uR1_perm,a2,"outflow"),[3 2 1]); 

        uR2_perm = permute(uR2,[3 2 1]); 
        uR21 = permute(WENO_pages(uR2_perm,a1,"outflow"),[3 2 1]);
        uR22 = permute(WENO_pages(uR2_perm,a2,"outflow"),[3 2 1]); 


        pL1_perm = permute(pL1,[3 2 1]); 
        pL11 = permute(WENO_pages(pL1_perm,a1,"outflow"),[3 2 1]);
        pL12 = permute(WENO_pages(pL1_perm,a2,"outflow"),[3 2 1]); 

        pL2_perm = permute(pL2,[3 2 1]); 
        pL21 = permute(WENO_pages(pL2_perm,a1,"outflow"),[3 2 1]);
        pL22 = permute(WENO_pages(pL2_perm,a2,"outflow"),[3 2 1]); 

        pR1_perm = permute(pR1,[3 2 1]); 
        pR11 = permute(WENO_pages(pR1_perm,a1,"outflow"),[3 2 1]);
        pR12 = permute(WENO_pages(pR1_perm,a2,"outflow"),[3 2 1]); 

        pR2_perm = permute(pR2,[3 2 1]); 
        pR21 = permute(WENO_pages(pR2_perm,a1,"outflow"),[3 2 1]);
        pR22 = permute(WENO_pages(pR2_perm,a2,"outflow"),[3 2 1]); 


    % Flux on the left interface (outflow)
    rhoim1_1  = [rhoR11(1,:,:);rhoR11(1:Nx-1,:,:)];
    rhoi_1    = rhoL11;
    rhoim1_2  = [rhoR12(1,:,:);rhoR12(1:Nx-1,:,:)];
    rhoi_2    = rhoL12;
    rhoim1_3  = [rhoR21(1,:,:);rhoR21(1:Nx-1,:,:)];
    rhoi_3    = rhoL21;
    rhoim1_4  = [rhoR22(1,:,:);rhoR22(1:Nx-1,:,:)];
    rhoi_4    = rhoL22;

    uim1_1  = [uR11(1,:,:);uR11(1:Nx-1,:,:)];
    ui_1    = uL11;
    uim1_2  = [uR12(1,:,:);uR12(1:Nx-1,:,:)];
    ui_2    = uL12;
    uim1_3  = [uR21(1,:,:);uR21(1:Nx-1,:,:)];
    ui_3    = uL21;
    uim1_4  = [uR22(1,:,:);uR22(1:Nx-1,:,:)];
    ui_4    = uL22;


    pim1_1  = [pR11(1,:,:);pR11(1:Nx-1,:,:)];
    pi_1    = pL11;
    pim1_2  = [pR12(1,:,:);pR12(1:Nx-1,:,:)];
    pi_2    = pL12;
    pim1_3  = [pR21(1,:,:);pR21(1:Nx-1,:,:)];
    pi_3    = pL21;
    pim1_4  = [pR22(1,:,:);pR22(1:Nx-1,:,:)];
    pi_4    = pL22;


    Eim1_1 = pim1_1 ./ (gamma- 1) + 0.5 * rhoim1_1 .* uim1_1.^2;
    Ei_1 = pi_1 ./ (gamma - 1) + 0.5 * rhoi_1 .* ui_1.^2;
    Eim1_2 = pim1_2 ./ (gamma - 1) + 0.5 * rhoim1_2 .* uim1_2.^2;
    Ei_2 = pi_2 ./ (gamma - 1) + 0.5 * rhoi_2 .* ui_2.^2;
    Eim1_3 = pim1_3 ./ (gamma - 1) + 0.5 * rhoim1_3 .* uim1_3.^2;
    Ei_3 = pi_3 ./ (gamma - 1) + 0.5 * rhoi_3 .* ui_3.^2;
    Eim1_4 = pim1_4 ./ (gamma - 1) + 0.5 * rhoim1_4 .* uim1_4.^2;
    Ei_4 = pi_4 ./ (gamma - 1) + 0.5 * rhoi_4 .* ui_4.^2;


    lambda1 = max(abs(ui_1)+real(sqrt(gamma.*pi_1./rhoi_1)),...
        abs(uim1_1)+real(sqrt(gamma.*pim1_1./rhoim1_1)));

    lambda2 = max(abs(ui_2)+real(sqrt(gamma.*pi_2./rhoi_2)),...
        abs(uim1_2)+real(sqrt(gamma.*pim1_2./rhoim1_2)));

    lambda3= max(abs(ui_3)+real(sqrt(gamma.*pi_3./rhoi_3)),...
        abs(uim1_3)+real(sqrt(gamma.*pim1_3./rhoim1_3)));

    lambda4 = max(abs(ui_4)+real(sqrt(gamma.*pi_4./rhoi_4)),...
        abs(uim1_4)+real(sqrt(gamma.*pim1_4./rhoim1_4)));
 

    F_rho_im12_1 = (rhoi_1.*ui_1 + rhoim1_1.*uim1_1)/2 - 0.5*lambda1.*(rhoi_1-rhoim1_1);
    F_rho_im12_2 = (rhoi_2.*ui_2 + rhoim1_2.*uim1_2)/2 - 0.5*lambda2.*(rhoi_2-rhoim1_2);
    F_rho_im12_3 = (rhoi_3.*ui_3 + rhoim1_3.*uim1_3)/2 - 0.5*lambda3.*(rhoi_3-rhoim1_3);
    F_rho_im12_4 = (rhoi_4.*ui_4 + rhoim1_4.*uim1_4)/2 - 0.5*lambda4.*(rhoi_4-rhoim1_4);
    
    F_rhou_im12_1 = (rhoi_1.*ui_1.^2+pi_1 + rhoim1_1.*uim1_1.^2+pim1_1)/2 - 0.5*lambda1.*(rhoi_1.*ui_1-rhoim1_1.*uim1_1);
    F_rhou_im12_2 = (rhoi_2.*ui_2.^2+pi_2 + rhoim1_2.*uim1_2.^2+pim1_2)/2 - 0.5*lambda2.*(rhoi_2.*ui_2-rhoim1_2.*uim1_2);
    F_rhou_im12_3 = (rhoi_3.*ui_3.^2+pi_3 + rhoim1_3.*uim1_3.^2+pim1_3)/2 - 0.5*lambda3.*(rhoi_3.*ui_3-rhoim1_3.*uim1_3);
    F_rhou_im12_4 = (rhoi_4.*ui_4.^2+pi_4 + rhoim1_4.*uim1_4.^2+pim1_4)/2 - 0.5*lambda4.*(rhoi_4.*ui_4-rhoim1_4.*uim1_4);


    F_E_im12_1 = (ui_1.*(Ei_1 + pi_1) + uim1_1.*(Eim1_1 + pim1_1))/2 - 0.5*lambda1.*(Ei_1-Eim1_1);
    F_E_im12_2 = (ui_2.*(Ei_2 + pi_2) + uim1_2.*(Eim1_2 + pim1_2))/2 - 0.5*lambda2.*(Ei_2-Eim1_2);
    F_E_im12_3 = (ui_3.*(Ei_3 + pi_3) + uim1_3.*(Eim1_3 + pim1_3))/2 - 0.5*lambda3.*(Ei_3-Eim1_3);
    F_E_im12_4 = (ui_2.*(Ei_4 + pi_4) + uim1_4.*(Eim1_4 + pim1_4))/2 - 0.5*lambda4.*(Ei_4-Eim1_4);



    Fim12_1 = [F_rho_im12_1; F_rhou_im12_1; F_E_im12_1];
    Fim12_2 = [F_rho_im12_2; F_rhou_im12_2; F_E_im12_2];
    Fim12_3 = [F_rho_im12_3; F_rhou_im12_3; F_E_im12_3];
    Fim12_4 = [F_rho_im12_4; F_rhou_im12_4; F_E_im12_4];



    % Flux on the right interface (outflow)
    rhoip1_1  = [rhoL11(2:Nx,:,:);rhoL11(Nx,:,:)];
    rhoi_1    = rhoR11;
    rhoip1_2  = [rhoL12(2:Nx,:,:);rhoL12(Nx,:,:)];
    rhoi_2    = rhoR12;
    rhoip1_3  = [rhoL21(2:Nx,:,:);rhoL21(Nx,:,:)];
    rhoi_3    = rhoR21;
    rhoip1_4  = [rhoL22(2:Nx,:,:);rhoL22(Nx,:,:)];
    rhoi_4    = rhoR22;

    uip1_1  = [uL11(2:Nx,:,:);uL11(Nx,:,:)];
    ui_1    = uR11;
    uip1_2  = [uL12(2:Nx,:,:);uL12(Nx,:,:)];
    ui_2    = uR12;
    uip1_3  = [uL21(2:Nx,:,:);uL21(Nx,:,:)];
    ui_3    = uR21;
    uip1_4  = [uL22(2:Nx,:,:);uL22(Nx,:,:)];
    ui_4    = uR22;


    pip1_1  = [pL11(2:Nx,:,:);pL11(Nx,:,:)];
    pi_1    = pR11;
    pip1_2  = [pL12(2:Nx,:,:);pL12(Nx,:,:)];
    pi_2    = pR12;
    pip1_3  = [pL21(2:Nx,:,:);pL21(Nx,:,:)];
    pi_3    = pR21;
    pip1_4  = [pL22(2:Nx,:,:);pL22(Nx,:,:)];
    pi_4    = pR22;


    Eip1_1 = pip1_1 ./ (gamma - 1) + 0.5 * rhoip1_1 .* uip1_1.^2;
    Ei_1 = pi_1 ./ (gamma - 1) + 0.5 * rhoi_1 .* ui_1.^2;
    Eip1_2 = pip1_2 ./ (gamma - 1) + 0.5 * rhoip1_2 .* uip1_2.^2;
    Ei_2 = pi_2 ./ (gamma - 1) + 0.5 * rhoi_2 .* ui_2.^2;
    Eip1_3 = pip1_3 ./ (gamma - 1) + 0.5 * rhoip1_3 .* uip1_3.^2;
    Ei_3 = pi_3 ./ (gamma - 1) + 0.5 * rhoi_3 .* ui_3.^2;
    Eip1_4 = pip1_4 ./ (gamma - 1) + 0.5 * rhoip1_4.* uip1_4.^2;
    Ei_4 = pi_4 ./ (gamma - 1) + 0.5 * rhoi_4 .* ui_4.^2;



    lambda1 = max(abs(ui_1)+real(sqrt(gamma.*pi_1./rhoi_1)),...
        abs(uip1_1)+real(sqrt(gamma.*pip1_1./rhoip1_1)));

    lambda2 = max(abs(ui_2)+real(sqrt(gamma.*pi_2./rhoi_2)),...
        abs(uip1_2)+real(sqrt(gamma.*pip1_2./rhoip1_2)));

    lambda3 = max(abs(ui_3)+real(sqrt(gamma.*pi_3./rhoi_3)),...
        abs(uip1_3)+real(sqrt(gamma.*pip1_3./rhoip1_3)));

    lambda4 = max(abs(ui_4)+real(sqrt(gamma.*pi_4./rhoi_4)),...
        abs(uip1_4)+real(sqrt(gamma.*pip1_4./rhoip1_4)));



    F_rho_ip12_1 = (rhoi_1.*ui_1 + rhoip1_1.*uip1_1)/2 - 0.5*lambda1.*(rhoip1_1-rhoi_1);
    F_rho_ip12_2 = (rhoi_2.*ui_2 + rhoip1_2.*uip1_2)/2 - 0.5*lambda2.*(rhoip1_2-rhoi_2);
    F_rho_ip12_3 = (rhoi_3.*ui_3 + rhoip1_3.*uip1_3)/2 - 0.5*lambda3.*(rhoip1_3-rhoi_3);
    F_rho_ip12_4 = (rhoi_4.*ui_4 + rhoip1_4.*uip1_4)/2 - 0.5*lambda4.*(rhoip1_4-rhoi_4);
    
    F_rhou_ip12_1 = (rhoi_1.*ui_1.^2+pi_1 + rhoip1_1.*uip1_1.^2+pip1_1)/2 - 0.5*lambda1.*(rhoip1_1.*uip1_1-rhoi_1.*ui_1);
    F_rhou_ip12_2 = (rhoi_2.*ui_2.^2+pi_2 + rhoip1_2.*uip1_2.^2+pip1_2)/2 - 0.5*lambda2.*(rhoip1_2.*uip1_2-rhoi_2.*ui_2);
    F_rhou_ip12_3 = (rhoi_3.*ui_3.^2+pi_3 + rhoip1_3.*uip1_3.^2+pip1_3)/2 - 0.5*lambda3.*(rhoip1_3.*uip1_3-rhoi_3.*ui_3);
    F_rhou_ip12_4 = (rhoi_4.*ui_4.^2+pi_4 + rhoip1_4.*uip1_4.^2+pip1_4)/2 - 0.5*lambda4.*(rhoip1_4.*uip1_4-rhoi_4.*ui_4);


    F_E_ip12_1 = (ui_1.*(Ei_1 + pi_1) + uip1_1.*(Eip1_1 + pip1_1))/2 - 0.5*lambda1.*(Eip1_1-Ei_1);
    F_E_ip12_2 = (ui_2.*(Ei_2 + pi_2) + uip1_2.*(Eip1_2 + pip1_2))/2 - 0.5*lambda2.*(Eip1_2-Ei_2);
    F_E_ip12_3 = (ui_3.*(Ei_3 + pi_3) + uip1_3.*(Eip1_3 + pip1_3))/2 - 0.5*lambda3.*(Eip1_3-Ei_3);
    F_E_ip12_4 = (ui_2.*(Ei_4 + pi_4) + uip1_4.*(Eip1_4 + pip1_4))/2 - 0.5*lambda4.*(Eip1_4-Ei_4);



    Fip12_1 = [F_rho_ip12_1; F_rhou_ip12_1; F_E_ip12_1];
    Fip12_2 = [F_rho_ip12_2; F_rhou_ip12_2; F_E_ip12_2];
    Fip12_3 = [F_rho_ip12_3; F_rhou_ip12_3; F_E_ip12_3];
    Fip12_4 = [F_rho_ip12_4; F_rhou_ip12_4; F_E_ip12_4];

    dUdt = reshape(-1/dx*((Fip12_1+Fip12_2+Fip12_3+Fip12_4)/4-...
        (Fim12_1+Fim12_2+Fim12_3+Fim12_4)/4), [],1);

end