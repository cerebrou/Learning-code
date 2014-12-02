clear; clc;
initpath;

learningType = 'dw';
%%
nP = 10;
nTrain = 50;
nTest = 100;
nOut = 0;
varmax = 0.05;
dFeat = 31;

%%
rng(0,'twister');
source = makeSourcePoints(nP, varmax);
trainPoints = makePoints(source, nTrain, 0);
% showPoints(source.xy);
% showPointsBoth(source, trainPoints{1});
refGraph = pointsToGraph(source);

%% Learn
[graphs, matches] = makeGraphs(trainPoints);

%%
nPair = nTrain;
nMatch = nP^2;
[L12(:,1), L12(:,2)] = find(ones(nMatch,nMatch));

for iPair = 1 : nPair
    M = makeAff(graphs{1}, graphs{iPair},[]);
    M(~triu(ones(nP^2))) = 0;
    vec{iPair} = M(M~=0);    
    [i,j] = find(M~=0);
    spc_i_indices{iPair} = i;
    spc_j_indices{iPair} = j;
    labels{iPair} = kron((1:nP)', ones(nP,1));
    nodes{iPair} = repmat((1:nP)',[nP,1]);
    nF(iPair) = nP;
end

learnSpectralMatching_Cars(vec, spc_i_indices, spc_j_indices, labels, nodes, nF);
% learnSpectralMatching_dw(vec, spc_i_indices, spc_j_indices, labels, nodes, nF);