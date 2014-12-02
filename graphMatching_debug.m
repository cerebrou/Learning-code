function [y,score] = graphMatching_debug(g1, g2,ytrue)
% INPUT
% g1, g2
% g1.fV: (dV by nV) matrix
% g1.fE: (dE by nE) matrix
% g1.edges: (2 by nE) matrix
% g1.nNodes
%
% OUTPUT
% y

% affinityMatrix = g1.fE' * g2.fE;
% affinityMatrix(logical(eye(size(affinityMatrix)))) = g1.fV .* g2.fV;
% yraw = RRWM(affinityMatrix);
% y = hungarian(yram);

nP = g1.nNodes;
nQ = g2.nNodes;
affinityMatrix = zeros(nP*nQ,nP*nQ);
for iEdgeModel = 1 : size(g1.edges,1)
    for iEdgeTarget = 1 : size(g2.edges,1)
        i = g1.edges(iEdgeModel,1);
        j = g1.edges(iEdgeModel,2);
        a = g2.edges(iEdgeTarget,1);
        b = g2.edges(iEdgeTarget,2);
        affinityMatrix(i+(a-1)*nP,j+(b-1)*nP) = ...
            g1.fE(:,iEdgeModel)' * g2.fE(:,iEdgeTarget);
    end
end
affOrig = affinityMatrix;
affinityMatrix = affinityMatrix - min(affinityMatrix(:));

%
E12 = ones(nP,nQ);
[L12(:,1) L12(:,2)] = find(E12);
[group1 group2] = make_group12(L12);
%
yraw = RRWM(double(affinityMatrix), group1, group2);
ysol = greedyMapping(yraw, group1, group2);
score = ysol(:)'*affOrig*ysol(:)
augScore = ysol(:)'*affinityMatrix*ysol(:)
ysol = reshape(ysol, [nP,nP]);
[y, ~] = find(ysol');

ytruebin = zeros(nP*nQ,1);
ytruebin( (1:nP) + ((ytrue)-1)*nP ) = 1;
trueScore = ytruebin(:)'*affOrig*ytruebin(:)
trueAugScore = ytruebin(:)'*affinityMatrix*ytruebin(:)
end