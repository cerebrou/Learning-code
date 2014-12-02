for iExp = 1 : 30
    iExp
    [~,i] = sort(accuracy(iExp,:));
    figure(1), clf, plot(accuracy(iExp,i));
    figure(2), clf, plot(score(iExp,i) / max(score(iExp,:)),'b');
    hold on; plot(scoresparseS(iExp,i) / max(scoresparseS(iExp,:)) ,'r'); hold off;
    pause;
end