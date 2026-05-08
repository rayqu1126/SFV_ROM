function [y_ids, z_ids, q12_y_ids, q34_y_ids, q13_z_ids, q24_z_ids, map_q] = HR_ids_2s(ids,Ny)

% Step 0: decode structure
cell_ids = ceil(ids/4);
quad_id  = mod(ids-1,4) + 1;

[i, j] = ind2sub([Ny, Ny], cell_ids);

% periodic neighbors
iL = mod(i-2, Ny) + 1;
iR = mod(i, Ny) + 1;

jD = mod(j-2, Ny) + 1;
jU = mod(j, Ny) + 1;

% ---- 9-point tensor stencil per cell ----
i_all = [
    iL; i; iR;   % x-left, center, right for each row
    iL; i; iR;
    iL; i; iR
];

j_all = [
    jD; jD; jD;
    j; j; j;
    jU; jU; jU
];

y_ids = unique(i_all);
z_ids = unique(j_all);
% 
mask_q12 = (quad_id == 1 | quad_id == 2);
mask_q34 = (quad_id == 3 | quad_id == 4);

y_q12 = unique(i(mask_q12));
y_q34 = unique(i(mask_q34));
[~, q12_y_ids] = ismember(y_q12, y_ids);
[~, q34_y_ids] = ismember(y_q34, y_ids);


mask_q13 = (quad_id == 1 | quad_id == 3);
mask_q24 = (quad_id == 2 | quad_id == 4);

z_q13 = unique(j(mask_q13));
z_q24 = unique(j(mask_q24));
[~, q13_z_ids] = ismember(z_q13, z_ids);
[~, q24_z_ids] = ismember(z_q24, z_ids);

lin = sub2ind([Ny,Ny], i, j);

% map_q(cell_index, quad_type) = position in ids
map_q = zeros(Ny^2, 4);

for k = 1:length(ids)
    map_q(lin(k), quad_id(k)) = k;
end

end
