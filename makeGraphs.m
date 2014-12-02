function [graphs, matches] = makeGraphs(points)

nSets = length(points);
graphs = cell(nSets,1);
matches = cell(nSets,1);
for iSet = 1 : nSets
    graphs{iSet} = pointsToGraph(points{iSet});
    matches{iSet} = points{iSet}.match(:);
end

end