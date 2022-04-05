%% This code is for validation of expected birth and death values

clc;
close all;
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Generate one network
p = 10;                              % number of nodes
q = p*(p-1)/2;                       % number of edges
m0 = p-1;                            % number of connected components - 1
m1 = (p-1)*(p-2)/2;                  % number of cycles
G = betarnd(1, 1, [1 q]);            % given network

n=15;                                % number of networks to be simulated
G1 = zeros(n,q);                     % initialization for simulated networks 
for i=1:n
    G1(i,:) = G + normrnd(0, 0.02, [1 q]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculate birth and death values of G directly using existing method

upper_tri_vec = G;
C = zeros(p,p);
C(logical(triu(ones(size(C)), 1))) = upper_tri_vec;
C = C + C.' + eye(p);            % calculate weighted adjacency matrix
    
b0 = conncomp_birth(C).';        % compute a set of increasing birth values
birth_val = b0(3,:);             % store original birth values

G_sort = sort(G);
death_val = G_sort(~ismember(G_sort,birth_val)); % store original death values

%% Empirical estimation of time-instances for group G1
b_G1 = zeros(n,m0);                  % initialize birth instances for group G1
d_G1 = zeros(n,m1);                  % initialize death instances for group G1
for i=1:n
    upper_tri_vec = G1(i,:);
    C = zeros(p,p);
    C(logical(triu(ones(size(C)), 1))) = upper_tri_vec;
    C = C + C.' + eye(p);            % calculate weighted adjacency matrix
    
    b0 = conncomp_birth(C).';        % compute a set of increasing birth values
    b_G1(i,:) = find(ismember(sort(upper_tri_vec),b0(3,:)));  % store birth instances
    d_G1(i,:) = find(~ismember(sort(upper_tri_vec),b0(3,:))); % store death instances
end

med_b_G1 = mean(b_G1);             % find mean of birth instances for group G1
med_d_G1 = mean(d_G1);             % find mean of death instances for group G1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Estimate birth and death values
[f,z] = ecdf(mean(G1));

% estimate birth values for first network x
r = (med_b_G1)/(q+1);
inv_F_b_x = zeros(1,m0);
for i=1:m0
    pos = find(r(i)<=f,1);
    inv_F_b_x(i) = z(pos);
end

birth_mean = zeros(1,m0);
birth_var = zeros(1,m0);
birth_low_CI = zeros(1,m0);
birth_up_CI = zeros(1,m0);
[kde,xi] = ksdensity(mean(G1), inv_F_b_x);

for i=1:m0
    birth_mean(i) = inv_F_b_x(i);
    birth_var(i) = (r(i) * (1-r(i)))/((q+1)*(kde(i))^2);
    birth_low_CI(i) = birth_mean(i) - 1.96*birth_var(i);
    birth_up_CI(i) = birth_mean(i) + 1.96*birth_var(i);
end

% estimate death values for first network x
r = (med_d_G1)/(q+1);
inv_F_d_x = zeros(1,m1);
for i=1:m1
    pos = find(r(i)<=f,1);
    inv_F_d_x(i) = z(pos);
end

death_mean = zeros(1,m1);
death_var = zeros(1,m1);
death_low_CI = zeros(1,m1);
death_up_CI = zeros(1,m1);
[kde,xi] = ksdensity(mean(G1), inv_F_d_x);

for i=1:m1
    death_mean(i) = inv_F_d_x(i);
    death_var(i) = (r(i) * (1-r(i)))/((q+1)*(kde(i))^2);
    death_low_CI(i) = death_mean(i) - 1.96*death_var(i);
    death_up_CI(i) = death_mean(i) + 1.96*death_var(i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot the ground truth vs estimated birth values

% eps = 0:0.0001:1;
% stairb = zeros(1,q);
% stairb(1) = 1;
% for i=2:length(eps)
%     if(logical(find(eps(i)==round(birth_val,4))))
%         stairb(i) = stairb(i-1) + 1;
%     else
%         stairb(i) = stairb(i-1);
%     end
% end
% 
% stairb_est(1) = 1;
% for i=2:length(eps)
%     if(logical(find(eps(i)==round(birth_mean,4))))
%         stairb_est(i) = stairb_est(i-1) + 1;
%     else
%         stairb_est(i) = stairb_est(i-1);
%     end
% end
% 
% stairb_est_low(1) = 1;
% for i=2:length(eps)
%     if(logical(find(eps(i)==round(birth_low_CI,4))))
%         stairb_est_low(i) = stairb_est_low(i-1) + 1;
%     else
%         stairb_est_low(i) = stairb_est_low(i-1);
%     end
% end
% 
% stairb_est_up(1) = 1;
% for i=2:length(eps)
%     if(logical(find(eps(i)==round(birth_up_CI,4))))
%         stairb_est_up(i) = stairb_est_up(i-1) + 1;
%     else
%         stairb_est_up(i) = stairb_est_up(i-1);
%     end
% end
% 
% stairs(eps,stairb, 'LineWidth',2, "color", "black")
% %xlim([0,q/(q+1)]);
% hold on
% stairs(eps,stairb_est, 'LineWidth',2, "color", "red")
% stairs(eps,stairb_est_low, 'LineWidth',2, "color", "blue")
% stairs(eps,stairb_est_up, 'LineWidth',2, "color", "blue")
% hold off
% set(gca,'FontSize',25);
% print(gcf,'birth_est','-dpng','-r300'); 

plot(birth_val,"color","black",'LineWidth',2,'LineStyle','-')
hold on
plot(birth_mean,"color","red",'LineWidth',2,'LineStyle','--')
plot(birth_low_CI,"color","blue",'LineWidth',2,'LineStyle','--')
plot(birth_up_CI,"color","blue",'LineWidth',2,'LineStyle','--')
hold off
set(gca,'FontSize',25);
print(gcf,'birth_est','-dpng','-r300'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot the ground truth vs estimated death values

plot(death_val,"color","black",'LineWidth',2,'LineStyle','-')
hold on
plot(death_mean,"color","red",'LineWidth',2,'LineStyle','--')
plot(death_low_CI,"color","blue",'LineWidth',2,'LineStyle','--')
plot(death_up_CI,"color","blue",'LineWidth',2,'LineStyle','--')
hold off 
set(gca,'FontSize',25);
print(gcf,'death_est','-dpng','-r300'); 
