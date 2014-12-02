function delta = lossHamming(param, y, ybar)
    delta = double(sum(y~=ybar)) / length(y);
end