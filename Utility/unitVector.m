function [unit_vec,norm_vec] = unitVector(vec)
    unit_vec = zeros(size(vec));
    norm_vec = (vec(1,:).^2 + vec(2,:).^2).^(0.5);
    nonzero_indx = find(norm_vec>0);
    unit_vec(:,nonzero_indx) = vec(:,nonzero_indx)./repmat(norm_vec(:,nonzero_indx),2,1);
end