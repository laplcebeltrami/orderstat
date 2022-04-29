% This code generates the validation table 2 for KS distance published in 
% Das, S., Anand, D.V., Chung, M.K. 2022 Topological data analysis for human 
% brain networks through order statistics arXiv:2204.02527.
%
% The code is downloaded from https://github.com/laplcebeltrami/orderstat

% (C) 2021 Soumya Das
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
n=6;                                 % number of networks in each group

tic
nrepeats=10;                         % number of iterations permutation test conducted
for iter=1:nrepeats
    G1 = betarnd(1, 1, [n q]);       % sample n networks for Group 1
    G2 = betarnd(1, 1, [n q]);       % sample n networks for Group 2
    b_G1 = zeros(n,m0);                  % initialize birth values for first group
    d_G1 = zeros(n,m1);                  % initialize death values for first group
    b_G2 = zeros(n,m0);                  % initialize birth values for second group
    d_G2 = zeros(n,m1);                  % initialize death values for second group
    for i=1:n
        upper_tri_vec = G1(i,:);
        M1 = zeros(p,p);
        M1(logical(triu(ones(size(M1)), 1))) = upper_tri_vec;
        M1 = M1 + M1.';            % calculate weighted adjacency matrix
        b0 = conncomp_birth(M1).';        % compute a set of increasing birth values
        b_G1(i,:) = find(ismember(sort(upper_tri_vec),b0(3,:)));  
        d_G1(i,:) = find(~ismember(sort(upper_tri_vec),b0(3,:))); 
       
        upper_tri_vec = G2(i,:);
        M2 = zeros(p,p);
        M2(logical(triu(ones(size(M2)), 1))) = upper_tri_vec;
        M2 = M2 + M2.';            % calculate weighted adjacency matrix
        b0 = conncomp_birth(M2).';        % compute a set of increasing birth values
        b_G2(i,:) = find(ismember(sort(upper_tri_vec),b0(3,:)));  
        d_G2(i,:) = find(~ismember(sort(upper_tri_vec),b0(3,:)));
    end

    % Observed statistics
    observed_distance = max(mean(b_G1)-mean(b_G2));  % computes KS-distance
    
    per_s = 100000;
    nG=n; % sample size in group I
    nP=n; % sample size in group II
    z=[b_G1; b_G2]; %[d_G1 d_G2];
    disp(iter)
    for i=1:per_s % each iteration gives a permutation
        store = randperm(nG+nP); % random permutation of data z.
        zper=z(store,:);
        xper=zper(1:nG,:);yper=zper((nG+1):(nG+nP),:); % permuted data is split into group 1 and 2  
        stat_s(i,iter) = max(mean(xper)-mean(yper));
     end
end
toc

% Online p-value computation
pvalues = online_pvalues(stat_s, observed_distance);
pvalend = pvalues(end)

