%% parameters
% nP = 10;
% nTest = 50;
% nOut = 0;
% varmax = 0.15;
% dFeat = 31;
methodid = {'SM','RRWM'};
nMethod = length(methodid);

%% generate points
% rng(0);
% source = makeSourcePoints(nP, varmax);
modelGraph = psiToGraph(model.w,nP);

% make problems
testPoints = makePoints(source, nTest, nOut);
vars = [0.001 0.01 0.1 1 10 100 1000 10000];
% vars = -pi/2 : pi/9 : pi/2;
nVar = length(vars);

% evaluate
acc = zeros(nMethod, nTest, nVar);
for iVar = 1 : nVar
    points = testPoints;
    disp(sprintf('%d-th variable',iVar));
    % rotate
%     scale= 1+rand();
    theta = 0;
    scale= vars(iVar)+rand();
%     theta = vars(iVar); %+rand();   % not the theta variance (~= [Oliver])
    Mrot = [cos(theta) -sin(theta) ; sin(theta) cos(theta) ];
    for iTest = 1 : nTest
        points{iTest}.xy = (Mrot*testPoints{iTest}.xy'*scale)';
    end
    
    for iTest = 1: nTest
        disp(sprintf('%d-th test',iTest));
        testGraph = pointsToGraph(points{iTest});
        % test to several methods
        for iMethod = 1 : nMethod
            y = graphMatching(modelGraph, testGraph, methodid{iMethod});
            acc(iMethod, iTest, iVar) = ...
                1 - sum(y~=points{iTest}.match(:))/model.nNode;
        end
    end
end

meanacc = squeeze(mean(acc,2));
% plot(vars,meanacc','linewidth',2);
semilogx(vars,meanacc','linewidth',2);
legend(methodid{1}, methodid{2});
xlabel('scale');
ylabel('accuracy');
title('graph matching');
ylim([0 1]);