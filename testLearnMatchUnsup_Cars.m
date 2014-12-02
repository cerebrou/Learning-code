
rand('state',sum(100*clock))

for nExp = 1 : 2 %lim(1):lim(2)

    %train_seq = randperm(100);
    
    train_seq = 2+randperm(30);
    
    train_seq = train_seq(1:10);
        
%     learnSpectralMatching_Cars;
    learnSpectralMatching_dw;

end