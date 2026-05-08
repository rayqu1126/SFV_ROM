function [y_ids,q1_ids,q2_ids,type_merged,perm] = HR_ids_1s(ids,Ny)
y_ids_temp = unique(ceil(ids/2));
y_ids = [];
for i = 1:length(y_ids_temp)
    stc2 = y_ids_temp(i);
    if stc2 == 1
        stc1 = Ny;
    else 
        stc1 = stc2-1;
    end

    if stc2 == Ny
        stc3 = 1;
    else
        stc3 = stc2+1;
    end
    y_ids = unique([y_ids;stc1;stc2;stc3]);
end

q1_y_ids = [];
q2_y_ids = [];
for i = 1:length(ids)
    if mod(ids(i),2) == 0
        q2_y_ids = [q2_y_ids;floor(ids(i)/2)];
    else
        q1_y_ids = [q1_y_ids;ceil(ids(i)/2)];
    end
end

q1_y_ids = unique(q1_y_ids);
q2_y_ids = unique(q2_y_ids);

[~, q1_ids] = ismember(q1_y_ids, y_ids);
[~, q2_ids] = ismember(q2_y_ids, y_ids);

ids_all = [q1_ids(:); q2_ids(:)];
type = [ones(length(q1_ids),1); 2*ones(length(q2_ids),1)]; 
[~, perm] = sortrows([ids_all, type], [1 2]);
type_merged = type(perm);

[~, perm] = sort(ids);
invperm = zeros(size(perm));
invperm(perm) = 1:length(ids);
perm = invperm;
end

