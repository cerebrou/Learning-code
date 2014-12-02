function affinityMatrix = makeAff(gm, g,  y)

nP = gm.nNode;
nQ = g.nNode;
if(isempty(y))
    iTrueMatches = [];
else
    iTrueMatches = zeros(nP*nQ,1);
    iTrueMatches((y-1)*nP + (1:nP)') = 1;
    iTrueMatches = find(iTrueMatches);
end

nP = gm.nNode;
nQ = g.nNode;
% affinityMatrix = gm.fE' * g.fE;   % wrong! nP^2 by nQ^2
% affinityMatrix = zeros(nP*nQ,nP*nQ);
% for iEdgeModel = 1 : size(gm.edges,1)
%     for iEdgeTarget = 1 : size(g.edges,1)
%         i = gm.edges(iEdgeModel,1);
%         j = gm.edges(iEdgeModel,2);
%         a = g.edges(iEdgeTarget,1);
%         b = g.edges(iEdgeTarget,2);
%         affinityMatrix(i+(a-1)*nP,j+(b-1)*nP) = ...
%             gm.fE(:,iEdgeModel)' * g.fE(:,iEdgeTarget);
%     end
% end
affinityMatrix = mexGetAff(double(gm.edges), double(g.edges), double(gm.fE), double(g.fE), nP, nQ);
affinityMatrix = affinityMatrix - min(affinityMatrix(:));
affinityMatrix((1:nP*nQ) + ((1:nP*nQ)-1)*nP*nQ) = ...
    affinityMatrix((1:nP*nQ) + ((1:nP*nQ)-1)*nP*nQ) + 1/nP;
affinityMatrix(iTrueMatches + (iTrueMatches-1)*nP*nQ) = ...
    affinityMatrix(iTrueMatches + (iTrueMatches-1)*nP*nQ) - 1/nP;
affinityMatrix = affinityMatrix - min(affinityMatrix(:));

end