%% WENO reconstruction for 3D array (page by page and row by row)
function U_tilde = WENO_pages(U,a,bc)
    [~, ~, Nz] = size(U);
    U_tilde = zeros(size(U));
    for k = 1:Nz
        U_tilde(:,:,k) = WENO_2Darray(U(:,:,k),a,bc);
    end
end