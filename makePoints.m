function points = makePoints(source, nSets, nOut)

nP = size(source.xy,1);
points = cell(nSets,1);
for iSet = 1 : nSets
%     id = randperm(nP+nOut);
    id = 1:nP;
    id = id(1:nP);
    points{iSet}.xy = rand(nP+nOut,2);
    points{iSet}.xy(id,:) = source.xy + bsxfun(@times, source.vars, randn(nP,2));
    points{iSet}.match = id;
end

end