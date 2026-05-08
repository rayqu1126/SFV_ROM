%% WENO reconstruction for 1D array
function U_tilde = WENO(U,a,bc)
U_tilde = zeros(size(U));
d0 = 0.5*(a*a + a - 1/12)/a;
d1 = 1 - d0;
Nx = length(U);
    for i = 1:Nx
        if i == 1
            if bc == "outflow"
                Ust1 = U(i);
            elseif bc == "periodic"
                Ust1 = U(Nx);
            else 
                error("Check BC");
            end
            Ust2 = U(i);
            Ust3 = U(i+1);
        elseif i == Nx
            Ust1 = U(i-1);
            Ust2 = U(i);
            if bc == "outflow"
                Ust3 = U(i);
            elseif bc == "periodic"
                Ust3 = U(1);
            end
        else
            Ust1 = U(i-1);
            Ust2 = U(i);
            Ust3 = U(i+1);
        end
     
        % smoothness indicator
        is0 = (Ust3 - Ust2)^2;
        is1 = (Ust2 - Ust1)^2;
        
        a0 = d0 / (eps + is0)^2;
        a1 = d1 / (eps + is1)^2;
        w0 = a0 / (a0 + a1);
        w1 = a1 / (a0 + a1);
        p0 = Ust2 + a*(Ust3 - Ust2);
        p1 = Ust2 + a*(Ust2 - Ust1);

        U_tilde(i) = w0 * p0 + w1 * p1;
    end
end
