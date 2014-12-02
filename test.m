nExp = 10;
for iExp = 1 : nExp

initpath;

learningType = 'dw-marius';
%%
nP = 10;
nTrain = 50;
nTest = 100;
nOut = 5;
varmax = 0.15;
dFeat = 31;

%%
% rng(0,'twister');
source = makeSourcePoints(nP, varmax);
trainPoints = makePoints(source, nTrain, 0);
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
testPoints = makePoints(source, nTest, nOut);

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
%         a(iExp) = mean(accTest);
    case 'dw-marius'
        wmatrix = zeros(nP,nP);
        wmatrix(~tril(ones(nP))) = model.w;
        wmatrix = wmatrix + wmatrix';
        for iTrain = 1 : nTrain
            trainGraph = pointsToGraph(trainPoints{iTrain});
            trainMatches = trainPoints{iTrain}.match;
            modelGraph = refGraph;
            modelGraph.fE = modelGraph.fE .* repmat(wmatrix(:)',dFeat,1);
            y = graphMatching(modelGraph, trainGraph, 'RRWM');
        %     [y,score] = graphMatching_debug(modelGraph, trainGraph, trainMatches);
            accTrain(iTrain) = 1 - sum(y~=trainPoints{iTrain}.match(:))/nP;
        end
        atrain(iExp) = mean(accTrain);
        for iTest = 1 : nTest
            testGraph = pointsToGraph(testPoints{iTest});
            testMatches = testPoints{iTest}.match;
            modelGraph = refGraph;
            modelGraph.fE = modelGraph.fE .* repmat(wmatrix(:)',dFeat,1);
            y = graphMatching(modelGraph, testGraph, 'RRWM');
        %     [y,score] = graphMatching_debug(modelGraph, testGraph, testMatches);
            accTest(iTest) = 1 - sum(y~=testPoints{iTest}.match(:))/nP;
        end
        atest(iExp) = mean(accTest);

end
% mean(a)
% hist(a)

end