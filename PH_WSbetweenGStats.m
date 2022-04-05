function betweenGStats = PH_WSbetweenGStats(lossMtx, nGroup1, nGroup2) 
%
% INPUT
% lossMtx          : matrix whose entries are pair-wise losses/distances
% nGroup1, nGroup2 : sample size


% OUTPUT
% observed : observed p-value

% between groups
% sum of pair-wise distances between groups
between = 0;
for i = 1:nGroup1
    for j = nGroup1+1:(nGroup1+nGroup2)
        between = between + lossMtx(i,j);
    end
end
denombetween = nGroup1*nGroup2;
betweenGStats = between/denombetween;