function visualizeEdges()

tag = 'house';
pFs = [1 100]; % frame index
nIn = [30 30] - 2; % randomly remove 2 nodes
parKnl = st('alg', 'cmum'); % type of affinity: only edge distance
[pars, algs] = gmPar(2);
wsSrc = cmumAsgSrc(tag, pFs, nIn, 'svL', 1);


F0s = cmumVdo(wsSrc, 'svL', 2);
Fs = F0s(wsSrc.pFs);
rows = 1; cols = 1;
Ax = iniAx(1, rows, cols, [400 * rows, 900 * cols], 'hGap', .1, 'wGap', .1);
parCor = st('cor', 'ln', 'mkSiz', 7, 'cls', {'y', 'b', 'g'});

box0s = imgBox(Fs);
[boxs, boxG, Lns, Sca] = boxCat('horz', box0s, 'gap', 0);
shImgBox(Fs, boxs);
shGph(Pt2s, Egs);

shAsgImg(Fs, gphs, asgFgmD, asgT, parCor, 'ax', Ax{1}, 'ord', 'n');
title('result of FGM-D');


end