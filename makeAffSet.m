function M = makeAffSet(graphs,refGraph,dFeat)

nGraph = length(graphs);
M = cell(nGraph, 1);
nP1 = refGraph.nNode; 
nP2 = graphs{1}.nNode;
nMatches = nP1*nP2;
for iGraph = 1 : nGraph
    M{iGraph} = zeros(nMatches, nMatches, dFeat);
    for iFeat = 1 : dFeat
        G1 = reshape(refGraph.fE(iFeat,:),[nP1,nP1]);
        G2 = reshape(graphs{iGraph}.fE(iFeat,:),[nP2,nP2]);
        M{iGraph}(:,:,iFeat) = repmat(G1,[nP2,nP2]) .* kron(G2,ones(nP1,nP1));
    end
end

end