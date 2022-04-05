clc;
close all;
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Generate two samples of networks
p = 10;                              % number of nodes
q = p*(p-1)/2;                       % number of edges
m0 = p-1;                            % number of connected components - 1
m1 = (p-1)*(p-2)/2;                  % number of cycles
n=12;                                % number of networks in each group
G1 = betarnd(5, 2, [n q]);           % sample n networks for Group 1
G2 = betarnd(5, 2, [n q]);           % sample n networks for Group 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Empirical estimation of time-instances for first group
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

med_b_G1 = median(b_G1);             % find median of birth instances for first group
med_d_G1 = median(d_G1);             % find median of death instances for first group

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Empirical estimation of time-instances for second group
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

med_b_G2 = median(b_G2);             % find median of birth instances for second group
med_d_G2 = median(d_G2);             % find median of death instances for second group

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Estimate birth and death values for first group
inv_F_b_x = zeros(n,m0);
inv_F_d_x = zeros(n,m1);
for j=1:n
    [f,z] = ecdf(G1(j,:));
    % estimate birth values
    r = (med_b_G1)/(q+1);
    for i=1:m0
        pos = find(r(i)<=f,1);
        inv_F_b_x(j,i) = z(pos);
    end

    % estimate death values
    r = (med_d_G1)/(q+1);
    for i=1:m1
        pos = find(r(i)<=f,1);
        inv_F_d_x(j,i) = z(pos);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimate birth and death values for second group
inv_F_b_y = zeros(n,m0);
inv_F_d_y = zeros(n,m1);
for j=1:n
    [f,z] = ecdf(G2(j,:));
    % estimate birth values
    r = (med_b_G2)/(q+1);
    for i=1:m0
        pos = find(r(i)<=f,1);
        inv_F_b_y(j,i) = z(pos);
    end

    % estimate death values
    r = (med_d_G2)/(q+1);
    for i=1:m1
        pos = find(r(i)<=f,1);
        inv_F_d_y(j,i) = z(pos);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Permutation test on x and y using stat
observed_distance = max(abs(mean(inv_F_b_x)-mean(inv_F_b_y))) + max(abs(mean(inv_F_d_x)-mean(inv_F_d_y))); % Observed Statistics
%sum((mean(inv_F_b_x)-mean(inv_F_b_y)).^2) + sum((mean(inv_F_d_x)-mean(inv_F_d_y)).^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nrepeats=10;
per_s = 100000;
nG=n; % sample size in group I
nP=n; % sample size in group II
z1=[inv_F_b_x; inv_F_b_y];
z2=[inv_F_d_x; inv_F_d_y];    % combine the data
tic
for j=1:nrepeats
    disp(j)
    for i=1:per_s %each iteration gives a permutation
        %random permutation of data z.
        store = randperm(nG+nP);
        z1per=z1(store,:); z2per=z2(store,:); 
        %permuted data is split into group 1 and 2  
        xper_b=z1per(1:nG,:);yper_b=z1per((nG+1):(nG+nP),:);
        xper_d=z2per(1:nG,:);yper_d=z2per((nG+1):(nG+nP),:);
        stat_s(i,j) = max(abs(mean(xper_b)-mean(yper_b))) + max(abs(mean(xper_d)-mean(yper_d)));
        % sum((mean(xper_b)-mean(yper_b)).^2) + sum((mean(xper_d)-mean(yper_d)).^2);
     end
end
toc
% % save stat_s1 stat_s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Online p-value computation
pvalues = online_pvalues(stat_s, observed_distance);
pvalend = pvalues(end)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot histogram to visualise distributions
%% FOLLOW THIS PLOT VISUALIZATION STYLE
figure;
subplot(1,2,1); 
histogram(stat_s,'FaceColor',[0.7 0.7 0.7],...
    'NumBins',10);
hold on
plot([observed_distance observed_distance],[0 10000],'--r','linewidth',2);
xlabel('Test Statistic')
set(gcf, 'Position', [400 400 600 250])
set(gca, 'fontsize',16)

subplot(1,2,2);
plot(pvalues,'k','linewidth',2);
xlim([0 100000]); ylim([0 1.0])
ylabel('p-values');
set(gcf, 'Position', [400 400 600 250])
set(gca, 'fontsize',16)
% print(gcf,'10_histlencycle_12.png','-dpng','-r300'); 


figure;
histogram(stat_s,'FaceColor',[0.7 0.7 0.7],...
    'NumBins',30);
hold on
plot([observed_distance observed_distance],[0 250000],'--r','linewidth',2);
xlabel('Test Statistic')
set(gcf, 'Position', [400 400 600 250])
set(gca, 'fontsize',16)
print(gcf,'maxgaphist2.png','-dpng','-r300'); 