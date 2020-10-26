function splitTestGroup (D, splitvar, splitval, refGrp, testvar)

% split data in D based on the split variable and value
nGrp = length(splitval);

for i=1:nGrp
    
    if i==refGrp
        continue;
    end;
    
    rangeGrp    = splitval{i};
    rangeRefGrp = splitval{refGrp};
    
    split = (D.(splitvar)>=rangeGrp(1)) & (D.(splitvar)<rangeGrp(2));
    refsplit = (D.(splitvar)>=rangeRefGrp(1)) & (D.(splitvar)<rangeRefGrp(2));

    D1 = getrow(D,split);
    D2 = getrow(D,refsplit);

    fprintf('Ref group vs group %2.2f-%2.2f\n',rangeGrp(1),rangeGrp(2));
    fprintf('--------------------------------\n\n');
    
    % parametric test
    ttest(D1.(testvar),D2.(testvar),2,'independent');

    % non-parametric test
    [p,h] = ranksum(D1.(testvar),D2.(testvar),'tail','both','alpha',0.05);
    fprintf('\nPVal: %1.3f, Hyp: %d\n\n',p,h);
end;
