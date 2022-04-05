function withinGStats = PH_WSwithinGStats(lossMtx, nGroup1, nGroup2) 

%
% INPUT
% lossMtx          : matrix whose entries are pair-wise losses/distances
% nGroup1, nGroup2 : sample size


% OUTPUT
% observed : observed p-value
%

% within groups
within = 0;
% sum of pair-wise distances within groups
for i = 1:nGroup1 % group 1
    for j = i + 1:nGroup1
        within = within + lossMtx(i, j);
    end
end

for i = nGroup1 + 1:nGroup1 + nGroup2 % group 2
    for j = i + 1:nGroup1 + nGroup2
        within = within + lossMtx(i, j);
    end
end
denomwithin = (nGroup1*(nGroup1-1) + nGroup2*(nGroup2-1))/2;
withinGStats = within/denomwithin;


