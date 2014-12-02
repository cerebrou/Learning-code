% k-fold cross validation

k = 10;
nValid = floor(nData/k);
nTrain = nData - nvalid;

orderAll = randperm(nTrain);
accVal = zeros(k,1);
for lambda = lamdaList
    for kk = 1 : k
        % Get training and valication set
        valSet = kfold(orderAll,k,kk);
        trainSet = orderAll(~ismember(orderAll,valSet));

        % Train
        

        % Validate
        accVal(kk) = 
    end
end