function model = learnGraph_l1(graphs, matches, nP, dFeat, refGraph, learningType)

parm.patterns = graphs;
parm.labels = matches;
parm.lossFn = @lossCB;
parm.nP = nP;
parm.dFeat = dFeat;
parm.e = 0.000000002;
parm.learningType = learningType;
parm.gm = refGraph;

switch(learningType)
    case 'non'
        model = [];
        return;
    case 'HARG'
        parm.dimension = dFeat * nP^2;
    case 'dw'
        parm.dimension = nP^2;
    case 'sw'
        parm.dimension = 2;
%     case 'sw-marius'
%         M = makeAffSet(graphs,refGraph,dFeat);
% %         M = makeAffSet(graphs,graphs{1});
%         model.w = trainWeight_marius(M,matches,nP,dFeat);
%         return;
end

% args = ' -c 1.0 -o 2 -v 3 -e 0.00002 -# 30';
% model = yumin_ssvm2(parm);
model = yumin_ssvm_cutting(parm);

end