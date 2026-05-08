%% WENO reconstruction for 3D array (page by page and row by row) at selected indices
function U_tilde = WENO_pages_HR(U,a,bc,ids)
    [~, N2, N3] = size(U);
    U_tilde = zeros(length(ids),N2,N3);
    for k = 1:N3
        U_tilde(:,:,k) = WENO_2Darray_HR(U(:,:,k),a,bc,ids);
    end
end