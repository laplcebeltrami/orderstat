% This code generates the validation table 2 for GH distance published in 
% Das, S., Anand, D.V., Chung, M.K. 2022 Topological data analysis for human 
% brain networks through order statistics arXiv:2204.02527.
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
n=6;                                 % number of networks in each group

tic
nrepeats=10;                         % number of iterations permutation test conducted
for iter=1:nrepeats
    G1 = betarnd(5, 2, [n q]);       % sample n networks for Group 1
    G2 = betarnd(1, 5, [n q]);       % sample n networks for Group 2
    C1 = cell(n, 1);
    C2 = cell(n, 1);
    for i=1:n
        upper_tri_vec = G1(i,:);
        M1 = zeros(p,p);
        M1(logical(triu(ones(size(M1)), 1))) = upper_tri_vec;
        M1 = M1 + M1.';            % calculate weighted adjacency matrix
        C1{i} = M1;
       
        upper_tri_vec = G2(i,:);
        M2 = zeros(p,p);
        M2(logical(triu(ones(size(M2)), 1))) = upper_tri_vec;
        M2 = M2 + M2.';            % calculate weighted adjacency matrix
        C2{i} = M2;
    end
    avg_C1 = zeros(p,p);
    avg_C2 = zeros(p,p);
    for i=1:n
        avg_C1 = avg_C1 + C1{i}/n;
        avg_C2 = avg_C2 + C2{i}/n;
    end
    
    SLA = PH_SLM(max(max(avg_C1)) - avg_C1);
    SLB = PH_SLM(max(max(avg_C2)) - avg_C2);
    
    % Observed statistics
    observed_distance = max(max(SLA-SLB));  % computes GH-distance
    
    per_s = 100000;
    nG=n; % sample size in group I
    nP=n; % sample size in group II
    z=[C1; C2];
    disp(iter)
    for i=1:per_s % each iteration gives a permutation
        store = randperm(nG+nP); % random permutation of data z.
        zper=z(store,:);
        xper=zper(1:nG,:);yper=zper((nG+1):(nG+nP),:); % permuted data is split into group 1 and 2  
        avg_xper = zeros(p,p);
        avg_yper = zeros(p,p);
        for k=1:n
            avg_xper = avg_xper + xper{k}/n;
            avg_yper = avg_yper + yper{k}/n;
        end
        SLA = PH_SLM(max(max(avg_xper)) - avg_xper);
        SLB = PH_SLM(max(max(avg_yper)) - avg_yper);
        stat_s(i,iter) = max(max(SLA-SLB));
     end
end
toc

% Online p-value computation
pvalues = online_pvalues(stat_s, observed_distance);
pvalend = pvalues(end)

