make2classProblem;

nTrain = n;
patterns = cell(nTrain,1);
y = cell(nTrain,1);
for iTrain = 1 : nTrain
    patterns{iTrain} = [X(iTrain,:), 1];
    y{iTrain} = labels(iTrain);
end
input.patterns = patterns;
input.labels = y;
input.dFeat = 3;

model = yumin_ssvm_cutting_test(input);
z = model.w;

xmat = cell2mat(input.patterns);
ymat = cell2mat(input.labels);
figure(1);  clf;
hh = {};
hh{1}=plot(xmat(ymat==1,1), xmat(ymat==1,2), 'o' );
hold all
hh{2}=plot(xmat(ymat==-1,1), xmat(ymat==-1,2), '*' );
legend('Class 1','Class 2');
legend([hh{1:2}],'Class 1','Class2');
m = z(1:2);
b = z(3);
grid = linspace(-.5,1,100);
plot(grid, -(b+m(1)*grid)/m(2));
drawnow;