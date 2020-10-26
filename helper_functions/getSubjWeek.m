function S=getSubjWeek (D,count)

%to exclude individuals with just one time point

[F,SubjN]=pivottable(D.SubjN,D.week,D.week,'length');
C=sum(~isnan(F),2);
goodSubj = SubjN(C>count);
S=getrow(D,ismember(D.SubjN,goodSubj));