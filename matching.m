function acc = matching(model, testPoints, methodid, y)


nP1 = model.nNode;
nP2 = testPoints.nNode;
t1 = model.t1;
t2 = [repmat(1:nP2,[1,nP2^2])', kron(repmat(1:nP2,[1,nP2]),ones(1,nP2))', kron(1:nP2,ones(1,nP2^2))' ];
t2 = t2(:,[3 2 1]);
nT = size(model.t1,2);

%generate features
feat1 = reshape(model.w, model.dFeat, nT);
feat2 = computeFeatureTheta(testPoints.xy, t2);

%number of nearest neighbors used for each triangle (results can be bad if
%too low)
nNN=nP1^2; %300;

%find the nearest neighbors
[inds, simm] = annquerySim(feat2, feat1, nNN);  % TODO: sim rather than dist?

%build the tensor
[i j k]=ind2sub([nP2,nP2,nP2],int32(inds));
tmp=repmat(1:nT,nNN,1);
indH = t1(:,tmp(:))' + [k(:)-1 j(:)-1 i(:)-1]*nP1;
valH = simm(:);
%initiatialize X
X=1/nP1*ones(nP1,nP2);
valH = valH - min(valH(:));

%power iteration
switch(methodid)
    case 'SM'
        [X2, score]=tensorMatching(X,[],[],[],[],indH,valH);
    case 'RRWM'
        X2=RRWHM(indH,valH,nP1,nP2);
end

%
E12 = ones(nP1,nP2);
[L12(:,1) L12(:,2)] = find(E12);
[group1 group2] = make_group12(L12);

%power iteration
yraw = X2;
ysol = greedyMapping(yraw(:)',group1,group2);
ysol = reshape(ysol, [nP1,nP2]);
[y,~] = find(ysol');
gsol = testPoints.xy(y,:);

acc = 1 - sum(y(:)~=testPoints.match(:))/nP1;
    
end