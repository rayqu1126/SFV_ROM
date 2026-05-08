%% WENO reconstruction for 2D array (row by row)
function U_tilde = WENO_2Darray(U,a,bc)
U_tilde = zeros(size(U));
d0 = 0.5*(a*a + a - 1/12)/a;
d1 = 1 - d0;
[N1, N2] = size(U);
for j = 1:N2
    for i = 1:N1
        if i == 1
            if bc == "outflow"
                Ust1 = U(i,j);
            elseif bc == "periodic"
                Ust1 = U(N1,j);
            else 
                error("Check BC");
            end
            Ust2 = U(i,j);
            Ust3 = U(i+1,j);
        elseif i == N1
            Ust1 = U(i-1,j);
            Ust2 = U(i,j);
            if bc == "outflow"
                Ust3 = U(i,j);
            elseif bc == "periodic"
                Ust3 = U(1,j);
            end
        else
            Ust1 = U(i-1,j);
            Ust2 = U(i,j);
            Ust3 = U(i+1,j);
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

        U_tilde(i,j) = w0 * p0 + w1 * p1;
    end
end
end
