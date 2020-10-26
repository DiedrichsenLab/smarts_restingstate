function splitTest (D, splitvar, splitval, testvar)

% split data in D based on the split variable and value
splitidx = D.(splitvar)<splitval;

D1 = getrow(D,splitidx);
D2 = getrow(D,~splitidx);

% parametric test
ttest(D1.(testvar),D2.(testvar),2,'independent');

% non-parametric test
[p,h] = ranksum(D1.(testvar),D2.(testvar),'tail','both','alpha',0.05);
fprintf('\n\nPVal: %1.3f, Hyp: %d\n',p,h);

keyboard;
