function [Wb Wd] = PH_decompose(W)
%[Wb, Wd] = PH_decompose(W)
%    
% The function computes the birth and death edge sets given in 
% Songdechakraiwut, T. Chung, M.K. 2020 Topological learning for brain 
% networks, arXiv: 2012.00675. If you are using any part of the code, 
% please reference the above paper.
%
%
% INPUT
% W : edge weight matrix
%
% OUTPUT
% Wb : birth edge set
% Wd : death edge set
%
%
% (C) 2020 Tananun Songdechakraiwut, Moo K. Chung
%          University of Wisconsin-Madison
%
%  Contact tananun@cs.wisc.edu or mkchung@wisc.edu
%  for support/permission with the codes 
%
% Update history
%     2020 November 11 created by Songdechakraiwut
%     2021 May 23 Modified Chung
%
%
%% Compute set of births and set of deaths

G1 = graph(W, 'upper', 'omitselfloops');

% birth edge set
birthMtx1 = conncomp_birth(W);
Wb=birthMtx1;

% death edge set
deathMtx1 = rmedge(G1, birthMtx1(:, 1), birthMtx1(:, 2)).Edges{:, :};
% sorting by weights in ascending order
deathMtx1 = sortrows(deathMtx1, 3);

Wd=deathMtx1;

