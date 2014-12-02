cbar = lines(nTrain);
range = 1 : nP^2*dFeat;
figure, hold on;
for iTrain = 1 : nTrain
    gperm = fitGraphToModel(graphs{iTrain},matches{iTrain});
    psi = graphToPsi(gperm);
    plot(psi(range),'color',cbar(iTrain,:));
end
plot(model.w(range)*100,'r');
hold off;


%%
figure, hold on;
gperm = fitGraphToModel(testGraph,testMatches);
psi = graphToPsi(gperm);
plot(psi,'color','r');
plot(model.w*100,'color','b');
hold off;

psi'*model.w