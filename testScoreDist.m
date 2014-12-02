% 
% nW = 100;
% dFeat = 31;
% wgroup = reshape(model.w, dFeat, 900);
% wgroup = sum(wgroup.^2);
% maxw = max(wgroup);
% % wlist = 0 : maxw/nW : maxw;
% wlist = maxw : -maxw/nW : 0;
% 
% model.w(model.w<0) = 0;
% modelGraph = phiToGraph(model.w,nP);
% scoreArch = zeros(nTest,nW);
% for iTest = 1 : nTest
%     testGraph = pointsToGraph(testPoints{iTest});
%     y = graphMatching(modelGraph, testGraph, 'RRWM');
%     matchedPoints.xy = testPoints{iTest}.xy(sort(y),:);
%     phi = graphToPhi(pointsToGraph(matchedPoints));
%     
%     scores = model.w .* phi(:);
%     totScore = sum(scores);
%     for iw = 1 : nW
%         biggroup = wgroup > wlist(iw);
%         bigel = kron(biggroup(:), ones(dFeat,1));
%         scoreArch(iTest,iw) = sum(scores(logical(bigel))) / totScore;
%     end
% end
% figure, plot(mean(scoreArch))
% figure, plot(scoreArch(1,:))


%%
dFeat = 31;
wgroup = reshape(model.w, dFeat, 900);
wgroup = sum(wgroup.^2);
maxw = max(wgroup);
[~,isort] = sort(wgroup,'descend');

modelGraph = phiToGraph(model.w,nP);
scoreArch = zeros(nTest,nW);
for iTest = 1 : nTest
    testGraph = pointsToGraph(testPoints{iTest});
    y = graphMatching(modelGraph, testGraph, 'RRWM');
    matchedPoints.xy = testPoints{iTest}.xy(sort(y),:);
    phi = graphToPhi(pointsToGraph(matchedPoints));
    
    scores = model.w .* phi(:);
%     scores = scores - min(scores(:));
    scoreGroup = sum(reshape(scores,dFeat,900));
    scoreGroupSort = scoreGroup(isort);
    totScore = sum(scores);
end
figure, plot(mean(scoreArch))
figure, plot(scoreArch(1,:))

