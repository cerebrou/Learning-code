% load('harg_wl21_92P3.mat')
load('harg_wl21_99P4.mat')
ww = reshape(model.w,31,900);
sumw = sum(ww.^2);
maxw = max(sumw);
big = sumw > maxw/30;
sum(big)
ww = ww .* repmat(big, 31, 1);
model.w = ww(:);
model.w(model.w<0) = 0;



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