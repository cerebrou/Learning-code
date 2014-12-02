initpath;

learningType = 'dw';
%%
nP = 10;
nTrain = 50;
nTest = 100;
nOut = 0;
varmax = 0.15;
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
model = learnGraph_l1(graphs, matches, nP, dFeat, refGraph, learningType);
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
        dimension = dFeat*nP^2;
        for iTrain = 1 : nTrain
            trainGraph = pointsToGraph(trainPoints{iTrain});
            trainMatches = trainPoints{iTrain}.match;
            modelGraph = phiToGraph(model.w(1:dimension),nP);
            y = graphMatching(modelGraph, trainGraph, 'RRWM');
        %     [y,score] = graphMatching_debug(modelGraph, trainGraph, trainMatches);
            accTrain(iTrain) = 1 - sum(y~=trainPoints{iTrain}.match(:))/nP;
        end
        mean(accTrain)
        for iTest = 1 : nTest
            testGraph = pointsToGraph(testPoints{iTest});
            testMatches = testPoints{iTest}.match;
            modelGraph = phiToGraph(model.w(1:dimension),nP);
            y = graphMatching(modelGraph, testGraph, 'RRWM');
        %     [y,score] = graphMatching_debug(modelGraph, testGraph, testMatches);
            accTest(iTest) = 1 - sum(y~=testPoints{iTest}.match(:))/nP;
        end
        mean(accTest)
    case 'dw'
        dimension = nP^2;
        for iTrain = 1 : nTrain
            trainGraph = pointsToGraph(trainPoints{iTrain});
            trainMatches = trainPoints{iTrain}.match;
            modelGraph = refGraph;
            modelGraph.fE = modelGraph.fE .* repmat(model.w(1:nP^2)',[dFeat,1]);
            y = graphMatching(modelGraph, trainGraph, 'RRWM');
        %     [y,score] = graphMatching_debug(modelGraph, trainGraph, trainMatches);
            accTrain(iTrain) = 1 - sum(y~=trainPoints{iTrain}.match(:))/nP;
        end
        mean(accTrain)
        for iTest = 1 : nTest
            testGraph = pointsToGraph(testPoints{iTest});
            testMatches = testPoints{iTest}.match;
            modelGraph = refGraph;
            modelGraph.fE = modelGraph.fE .* repmat(model.w(1:nP^2)',[dFeat,1]);
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

% figure, plot(sourcePsi(:))
% hold on; plot(model.w(:)*10,'r');
plot(model.w(1:dimension),'r');

ww = reshape(model.w(1:dimension),dFeat,nP^2);
figure, plot(sum(ww.^2,1));
on = kron( sum(ww.^2,1)>1e-7, ones(dFeat,1) );
model.w(~on) = 0;