randn('state',23432);
rand('state',3454);

n = 30;
d = 2; % 2D data is nice because it's easy to visualize

n1 = round(.5*n);   % number of -1 (this is data from one "class")
n2 = n - n1;        % number of +1 (this is data from the other "class")

% Generate data
mean1   = [-.5,1];
% mean2 = [.6,.2]; % can be separated
mean2   = [.1,.5]; % cannot be separated, so a bit more interesting
s       = .25;  % standard deviation
x1  = repmat(mean1,n1,1) + s*randn(n1,2);
x2  = repmat(mean2,n2,1) + s*randn(n2,2);

X = [x1; x2];
labels = [ ones(n1,1); -ones(n2,1) ];

% Plot the data:
figure(1); clf;
hh = {};
hh{1}=plot(x1(:,1), x1(:,2), 'o' );
hold all
hh{2}=plot(x2(:,1), x2(:,2), '*' );
legend('Class 1','Class 2');
xl = get(gca,'xlim');
yl = get(gca,'ylim');
legend([hh{1:2}],'Class 1','Class2');