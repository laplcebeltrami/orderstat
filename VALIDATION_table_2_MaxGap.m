% This code generates the validation table 2 for maximum gap statistic published 
% in Das, S., Anand, D.V., Chung, M.K. 2022 Topological data analysis for 
% human brain networks through order statistics arXiv:2204.02527.
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
n=6;                                % number of networks in each group

tic
nrepeats=10;                         % number of iterations permutation test conducted
for iter=1:nrepeats
    G1 = betarnd(1, 1, [n q]);           % sample n networks for Group 1
    G2 = betarnd(5, 2, [n q]);           % sample n networks for Group 2

    % Empirical estimation of time-instances for first group
    b_G1 = zeros(n,m0);                  % initialize birth instances for first group
    d_G1 = zeros(n,m1);                  % initialize death instances for first group
    for i=1:n
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
    b_G2 = zeros(n,m0);                  % initialize birth instances for second group
    d_G2 = zeros(n,m1);                  % initialize birth instances for second group
    for i=1:n
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
    inv_F_b_x = zeros(n,m0);
    inv_F_d_x = zeros(n,m1);
    for j=1:n
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
    inv_F_b_y = zeros(n,m0);
    inv_F_d_y = zeros(n,m1);
    for j=1:n
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

    % Observed statistics
    observed_distance = max(abs(mean(inv_F_b_x)-mean(inv_F_b_y))) + max(abs(mean(inv_F_d_x)-mean(inv_F_d_y)));

    per_s = 100000;
    nG=n; % sample size in group I
    nP=n; % sample size in group II
    z1=[inv_F_b_x; inv_F_b_y];
    z2=[inv_F_d_x; inv_F_d_y];    % combine the data
    disp(iter)
    for i=1:per_s %each iteration gives a permutation
        %random permutation of data z.
        store = randperm(nG+nP);
        z1per=z1(store,:); z2per=z2(store,:); 
        %permuted data is split into group 1 and 2  
        xper_b=z1per(1:nG,:);yper_b=z1per((nG+1):(nG+nP),:);
        xper_d=z2per(1:nG,:);yper_d=z2per((nG+1):(nG+nP),:);
        stat_s(i,iter) = max(abs(mean(xper_b)-mean(yper_b))) + max(abs(mean(xper_d)-mean(yper_d)));
     end
end
toc

% Online p-value computation
pvalues = online_pvalues(stat_s, observed_distance);
pvalend = pvalues(end)

% Plot histogram to visualise distributions
figure;
histogram(stat_s,'FaceColor',[0.7 0.7 0.7],...
    'NumBins',30);
hold on
plot([observed_distance observed_distance],[0 110000],'--r','linewidth',2);
xlabel('Test Statistic')
set(gcf, 'Position', [400 400 600 250])
set(gca, 'fontsize',16)
%print(gcf,'maxgaphist1.png','-dpng','-r300'); 
