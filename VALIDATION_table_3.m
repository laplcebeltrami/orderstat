% This code generates the validation table 3 for AUC (area under Betti curve) 
% statistic published in Das, S., Anand, D.V., Chung, M.K. 2022 Topological 
% data analysis for human brain networks through order statistics arXiv:2204.02527.
%
% The code is downloaded from https://github.com/laplcebeltrami/orderstat

% (C) 2021 Soumya Das, Moo K. Chung
%     University of Wisconsin-Madison
%
% Update history:
% Credited Nov 1, 2021  Das
% Edited   April 27, 2022 Chung

clc;
close all;
clear all;

%-----------
% Generate two samples of networks
p = 10;                              % number of nodes
q = p*(p-1)/2;                       % number of edges
m0 = p-1;                            % number of connected components - 1
m1 = (p-1)*(p-2)/2;                  % number of cycles
nn = 6;
nG1=nn;                              % number of networks in first group
nG2=nn;                              % number of networks in second group

G1 = betarnd(1, 1, [nG1 q]);           % sample n networks for Group 1
G2 = betarnd(5, 2, [nG2 q]);           % sample n networks for Group 2

% Empirical estimation of time-instances for first group
b_G1 = zeros(nG1,m0);                  % initialize birth instances for first group
d_G1 = zeros(nG1,m1);                  % initialize death instances for first group
for i=1:nG1
    upper_tri_vec = G1(i,:);
    C = zeros(p,p);
    C(logical(triu(ones(size(C)), 1))) = upper_tri_vec;
    C = C + C.' + eye(p);            % calculate weighted adjacency matrix
    
    b0 = conncomp_birth(C).';        % compute a set of increasing birth values
    b_G1(i,:) = find(ismember(sort(upper_tri_vec),b0(3,:)));  % store birth instances
    d_G1(i,:) = find(~ismember(sort(upper_tri_vec),b0(3,:))); % store death instances
end

mean_b_G1 = mean(b_G1);             % find mean of birth instances for first group
mean_d_G1 = mean(d_G1);             % find mean of death instances for first group

% Empirical estimation of time-instances for second group
b_G2 = zeros(nG2,m0);                  % initialize birth instances for second group
d_G2 = zeros(nG2,m1);                  % initialize birth instances for second group
for i=1:nG2
    upper_tri_vec = G2(i,:);
    C = zeros(p,p);
    C(logical(triu(ones(size(C)), 1))) = upper_tri_vec;
    C = C + C.' + eye(p);            % calculate weighted adjacency matrix
    
    b0 = conncomp_birth(C).';        % compute a set of increasing birth values
    b_G2(i,:) = find(ismember(sort(upper_tri_vec),b0(3,:)));  % store birth instances
    d_G2(i,:) = find(~ismember(sort(upper_tri_vec),b0(3,:))); % store death instances
end

mean_b_G2 = mean(b_G2);             % find mean of birth instances for second group
mean_d_G2 = mean(d_G2);             % find mean of death instances for second group

% Estimate birth and death values for first group
inv_F_b_x = zeros(nG1,m0);
inv_F_d_x = zeros(nG1,m1);
for j=1:nG1
    [f,z] = ecdf(G1(j,:));
    % Estimate birth values
    r = (mean_b_G1)/(q+1);
    for i=1:m0
        pos = find(r(i)<=f,1);
        inv_F_b_x(j,i) = z(pos);
    end

    % Estimate death values
    r = (mean_d_G1)/(q+1);
    for i=1:m1
        pos = find(r(i)<=f,1);
        inv_F_d_x(j,i) = z(pos);
    end
end

% Estimate birth and death values for second group
inv_F_b_y = zeros(nG2,m0);
inv_F_d_y = zeros(nG2,m1);
for j=1:nG2
    [f,z] = ecdf(G2(j,:));
     % Estimate birth values
     r = (mean_b_G2)/(q+1);
     for i=1:m0
         pos = find(r(i)<=f,1);
         inv_F_b_y(j,i) = z(pos);
     end

     % Estimate death values
     r = (mean_d_G2)/(q+1);
     for i=1:m1
         pos = find(r(i)<=f,1);
         inv_F_d_y(j,i) = z(pos);
     end
end

b1 = [zeros(size(inv_F_b_x,1),1) inv_F_b_x];
d1 = [zeros(size(inv_F_d_x,1),1) inv_F_d_x];
b2 = [zeros(size(inv_F_b_y,1),1) inv_F_b_y];
d2 = [zeros(size(inv_F_d_y,1),1) inv_F_d_y];
AUC1 = zeros(nG1,1);
for i=1:nG1
    for k=2:m0
        AUC1(i) = AUC1(i) + k*(b1(i,(k+1)) - b1(i,k));
    end
end
    
AUC2 = zeros(nG2,1);
for i=1:nG2
    for k=2:m0
        AUC2(i) = AUC2(i) + k*(b2(i,(k+1)) - b2(i,k));
    end
end

% Wilcoxon rank-sum test
[p,h] = ranksum(AUC1,AUC2);
p        % pvalue
h;       % decision
