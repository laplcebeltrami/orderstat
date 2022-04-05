function s = expected_topological_loss(x,y,med_b_G1,med_d_G1,med_b_G2,med_d_G2,q,m0,m1)
% function s = expected_topological_loss(x,y,med_b_G1,med_d_G1,med_b_G2,med_d_G2,q,m0,m1)
%
% This function computes the expected topological loss between two
% networks.
%
% Input: 
% x        : first network
% y        : second network
% med_b_G1 : empirical estimate of the birth instances of first network
% med_d_G1 : empirical estimate of the death instances of first network
% med_b_G2 : empirical estimate of the birth instances of second network
% med_d_G2 : empirical estimate of the death instances of second network
% q        : number of edges
% m0       : number of nodes - 1
% m1       : number of cycles
%
% Output:
% s        : Expected topological loss

% First network
[f1,z] = ecdf(x);

% estimate birth values for first network x
quantiles_b_x = (med_b_G1)/(q+1);
inv_F_b_x = zeros(1,m0);
for i=1:m0
    pos = find(quantiles_b_x(i)<=f1,1);
    inv_F_b_x(i) = z(pos);
end

% estimate death values for first network x
quantiles_d_x = (med_d_G1)/(q+1);
inv_F_d_x = zeros(1,m1);
for i=1:m1
    pos = find(quantiles_d_x(i)<=f1,1);
    inv_F_d_x(i) = z(pos);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Second network

[f2,z] = ecdf(y);

% estimate birth values for second network y
quantiles_b_y = (med_b_G2)/(q+1);
inv_F_b_y = zeros(1,m0);
for i=1:m0
    pos = find(quantiles_b_y(i)<=f2,1);
    inv_F_b_y(i) = z(pos);
end

% estimate death values for second network y
quantiles_d_y = (med_d_G2)/(q+1);
inv_F_d_y = zeros(1,m1);
for i=1:m1
    pos = find(quantiles_d_y(i)<=f2,1);
    inv_F_d_y(i) = z(pos);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s = sum((inv_F_b_x - inv_F_b_y) .^ 2) + sum((inv_F_d_x - inv_F_d_y) .^ 2);

%norm(inv_F_b_x - inv_F_b_y) + norm(inv_F_d_x-inv_F_d_y);

end