function g = pointsToGraph(points)
% INPUT
% points.xy: (nP by 2) matrix
% points.vars: (nP by 1) vector
%
% OUTPUT
% g.fV: (dV by nV) matrix
% g.fE: (dE by nE) matrix
% g.edges: (2 by nE) matrix
% g.nNodes

nP = size(points.xy,1);
g.fV = [];
% g.edges = combntns(1:nP,2);
g.edges = [repmat((1:nP)',nP,1), kron((1:nP)',ones(nP,1))];
nE = size(g.edges,1);
rawE = points.xy(g.edges(:,1),:) - points.xy(g.edges(:,2),:);
attrE = [sqrt(sum(rawE.^2,2)), atan2(rawE(:,2),rawE(:,1))]';
g.fE = makeFeatBin(attrE);
g.nNode = nP;