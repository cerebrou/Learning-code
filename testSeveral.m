nOut= 5;
for iExp = 1 : 30
    testPoints = makePoints(source, nTest, nOut);

        model.w = [0.5 0.5];
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
        a(iExp) = mean(accTest);
        
end