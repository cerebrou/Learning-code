function y = graphMatching(g1, g2, methodid)
% INPUT
% g1, g2
% g1.fV: (dV by nV) matrix
% g1.fE: (dE by nE) matrix
% g1.edges: (2 by nE) matrix
% g1.nNodes
%
% OUTPUT
% y

nP = g1.nNode;
nQ = g2.nNode;
% affinityMatrix = zeros(nP*nQ,nP*nQ);
% for iEdgeModel = 1 : size(g1.edges,1)
%     for iEdgeTarget = 1 : size(g2.edges,1)
%         i = g1.edges(iEdgeModel,1);
%         j = g1.edges(iEdgeModel,2);
%         a = g2.edges(iEdgeTarget,1);
%         b = g2.edges(iEdgeTarget,2);
%         affinityMatrix(i+(a-1)*nP,j+(b-1)*nP) = ...
%             g1.fE(:,iEdgeModel)' * g2.fE(:,iEdgeTarget);
%     end
% end
affinityMatrix = makeAff(g1,g2,[]);
affinityMatrix = affinityMatrix - min(affinityMatrix(:));
%
E12 = ones(nP,nQ);
[L12(:,1) L12(:,2)] = find(E12);
[group1 group2] = make_group12(L12);
%
switch(methodid)
    case 'SM'
        yraw = SM(double(affinityMatrix));
    case 'RRWM'
        yraw = RRWM(double(affinityMatrix), group1, group2);
end
% 
ysol = greedyMapping(yraw, group1, group2);
ysol = reshape(ysol, [nP,nQ]);
[y, ~] = find(ysol');

end