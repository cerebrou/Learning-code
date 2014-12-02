initpath;
global footpath
footpath = 'e:/codes/fgm-master';
addpath(genpath([footpath '/src']));
addpath(genpath([footpath '/lib']));

clear variables;
prSet(1);
learningType = 'HARG';

%% src parameter
nData = 111;
tag = 'house';
pFs = 1:nData; % frame index
nIn = 30*ones(nData,1); % randomly remove 2 nodes
parKnl = st('alg', 'cmum'); % type of affinity: only edge distance

%% algorithm parameter
[pars, algs] = gmPar(2);

%% src
rng(0,'twister');
CMUM = cmumHuman;
Pts = CMUM.(tag).XTs(pFs);
nTrain = round(nData/2);
nTest = nData - nTrain -1;
idx = randperm(nData);
idxTrain = idx(2:(nTrain+1));
idxTest = idx((nTrain+2):end);
source.xy = Pts{idx(1)}';
for iTrain = 1 : nTrain
    id = idxTrain(iTrain);
    % randomly pick inlier nodes
    n0 = size(Pts{id},2);
    ord = randperm(n0);
    ord = ord(1 : nIn(id));
    trainPoints{iTrain}.xy(ord,:) = Pts{id}';
    trainPoints{iTrain}.match = ord;
end
for iTest = 1 : nTest
    id = idxTest(iTest);
    % randomly pick inlier nodes
    n0 = size(Pts{id},2);
    ord = randperm(n0);
    ord = ord(1 : nIn(id));
    testPoints{iTest}.xy(ord,:) = Pts{id}';
    testPoints{iTest}.match = ord;
end

%%
nP = size(source.xy,1);
% nP = 10;
% nTrain = 50;
% nTest = 100;
nOut = 0;
% varmax = 0.05;
dFeat = 31;

%%
% source = makeSourcePoints(nP, varmax);
% trainPoints = makePoints(source, nTrain, 0);
% showPoints(source.xy);
% showPointsBoth(source, trainPoints{1});
refGraph = pointsToGraph(source);

%% Learn
[graphs, matches] = makeGraphs(trainPoints);
tic
model = learnGraph(graphs, matches, nP, dFeat,refGraph,learningType);
toc
model.nNode = nP;

sourceGraph = pointsToGraph(source);
sourcePsi = graphToPhi(sourceGraph);
%% Evaluate
accTrain = zeros(nTrain,1);
accTest = zeros(nTest,1);

% rotate
% scale= 0.5+rand();
% theta = 0.5*rand();
% Mrot = [cos(theta) -sin(theta) ; sin(theta) cos(theta) ];
% for iTest = 1 : nTest
%     testPoints{iTest}.xy = (Mrot*testPoints{iTest}.xy'*scale)';
% end

switch(learningType)
    case 'non'
        for iTrain = 1 : nTrain
            trainGraph = pointsToGraph(trainPoints{iTrain});
            trainMatches = trainPoints{iTrain}.match;
            y = graphMatching(refGraph, trainGraph, 'RRWM');
        %     [y,score] = graphMatching_debug(modelGraph, trainGraph, trainMatches);
            accTrain(iTrain) = 1 - sum(y~=trainPoints{iTrain}.match(:))/nP;
        end
        mean(accTrain)
        for iTest = 1 : nTest
            testGraph = pointsToGraph(testPoints{iTest});
            testMatches = testPoints{iTest}.match;
            y = graphMatching(refGraph, testGraph, 'RRWM');
        %     [y,score] = graphMatching_debug(modelGraph, testGraph, testMatches);
            accTest(iTest) = 1 - sum(y~=testPoints{iTest}.match(:))/nP;
        end
        mean(accTest)
    case 'HARG'
        for iTrain = 1 : nTrain
            trainGraph = pointsToGraph(trainPoints{iTrain});
            trainMatches = trainPoints{iTrain}.match;
            modelGraph = phiToGraph(model.w,nP);
            y = graphMatching(modelGraph, trainGraph, 'RRWM');
        %     [y,score] = graphMatching_debug(modelGraph, trainGraph, trainMatches);
            accTrain(iTrain) = 1 - sum(y~=trainPoints{iTrain}.match(:))/nP;
        end
        mean(accTrain)
        for iTest = 1 : nTest
            testGraph = pointsToGraph(testPoints{iTest});
            testMatches = testPoints{iTest}.match;
            modelGraph = phiToGraph(model.w,nP);
            y = graphMatching(modelGraph, testGraph, 'RRWM');
        %     [y,score] = graphMatching_debug(modelGraph, testGraph, testMatches);
            accTest(iTest) = 1 - sum(y~=testPoints{iTest}.match(:))/nP;
        end
        mean(accTest)
    case 'dw'
        for iTrain = 1 : nTrain
            trainGraph = pointsToGraph(trainPoints{iTrain});
            trainMatches = trainPoints{iTrain}.match;
            modelGraph = refGraph;
            modelGraph.fE = modelGraph.fE .* repmat(model.w(:)',[dFeat,1]);
            y = graphMatching(modelGraph, trainGraph, 'RRWM');
        %     [y,score] = graphMatching_debug(modelGraph, trainGraph, trainMatches);
            accTrain(iTrain) = 1 - sum(y~=trainPoints{iTrain}.match(:))/nP;
        end
        mean(accTrain)
        for iTest = 1 : nTest
            testGraph = pointsToGraph(testPoints{iTest});
            testMatches = testPoints{iTest}.match;
            modelGraph = refGraph;
            modelGraph.fE = modelGraph.fE .* repmat(model.w(:)',[dFeat,1]);
            y = graphMatching(modelGraph, testGraph, 'RRWM');
        %     [y,score] = graphMatching_debug(modelGraph, testGraph, testMatches);
            accTest(iTest) = 1 - sum(y~=testPoints{iTest}.match(:))/nP;
        end
        mean(accTest)
    case 'sw'
        for iTrain = 1 : nTrain
            trainGraph = pointsToGraph(trainPoints{iTrain});
            trainMatches = trainPoints{iTrain}.match;
            modelGraph = refGraph;
            modelGraph.fE(1:13,:) = modelGraph.fE(1:13,:) .* repmat(model.w(1),[13,nP^2]);
            modelGraph.fE(14:end,:) = modelGraph.fE(14:end,:) .* repmat(model.w(2),[18,nP^2]);
            y = graphMatching(modelGraph, trainGraph, 'RRWM');
        %     [y,score] = graphMatching_debug(modelGraph, trainGraph, trainMatches);
            accTrain(iTrain) = 1 - sum(y~=trainPoints{iTrain}.match(:))/nP;
        end
        mean(accTrain)
        for iTest = 1 : nTest
            testGraph = pointsToGraph(testPoints{iTest});
            testMatches = testPoints{iTest}.match;
            modelGraph = refGraph;
            modelGraph.fE(1:13,:) = modelGraph.fE(1:13,:) .* repmat(model.w(1),[13,nP^2]);
            modelGraph.fE(14:end,:) = modelGraph.fE(14:end,:) .* repmat(model.w(2),[18,nP^2]);
            y = graphMatching(modelGraph, testGraph, 'RRWM');
        %     [y,score] = graphMatching_debug(modelGraph, testGraph, testMatches);
            accTest(iTest) = 1 - sum(y~=testPoints{iTest}.match(:))/nP;
        end
        mean(accTest)
    case 'sw-marius'
        for iTrain = 1 : nTrain
            trainGraph = pointsToGraph(trainPoints{iTrain});
            trainMatches = trainPoints{iTrain}.match;
            modelGraph = refGraph;
            modelGraph.fE(1:13,:) = modelGraph.fE(1:13,:) .* repmat(model.w(1),[13,nP^2]);
            modelGraph.fE(14:end,:) = modelGraph.fE(14:end,:) .* repmat(model.w(2),[18,nP^2]);
            y = graphMatching(modelGraph, trainGraph, 'RRWM');
        %     [y,score] = graphMatching_debug(modelGraph, trainGraph, trainMatches);
            accTrain(iTrain) = 1 - sum(y~=trainPoints{iTrain}.match(:))/nP;
        end
        mean(accTrain)
        for iTest = 1 : nTest
            testGraph = pointsToGraph(testPoints{iTest});
            testMatches = testPoints{iTest}.match;
            modelGraph = refGraph;
            modelGraph.fE(1:13,:) = modelGraph.fE(1:13,:) .* repmat(model.w(1),[13,nP^2]);
            modelGraph.fE(14:end,:) = modelGraph.fE(14:end,:) .* repmat(model.w(2),[18,nP^2]);
            y = graphMatching(modelGraph, testGraph, 'RRWM');
        %     [y,score] = graphMatching_debug(modelGraph, testGraph, testMatches);
            accTest(iTest) = 1 - sum(y~=testPoints{iTest}.match(:))/nP;
        end
        mean(accTest)
end