function gbar = fitGraphToModel(g,y)
% INPUT
% g.fV: (dV by nV) matrix
% g.fE: (dE by nE) matrix
% g.edges: (2 by nE) matrix
% g.nNodes
% y: (nP by 1) matrix
%
% OUTPUT
% gbar.fV: (dV by nV) matrix
% gbar.fE: (dE by nE) matrix
% gbar.edges: (2 by nE) matrix

% disp('in fitGraphToModel');

nP = length(y);
nQ = g.nNode;
if(~isempty(g.fV))
    gbar.fV = g.fV(y);
else
    gbar.fV = [];
end

index = reshape(1:nQ^2,nQ,nQ);
index = index(y,y);
gbar.fE = g.fE(:,index(:));
gbar.edges = [repmat((1:nP)',nP,1), kron((1:nP)',ones(nP,1))];
gbar.nNodes = nP;

end