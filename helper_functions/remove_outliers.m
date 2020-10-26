case 'remove_outliers' 
    
D = load ('alldat.mat');
D = getrow (D, D.isGood==1);

[F,R]=pivottable([D.SubjN,D.week,D.state, D.task],[],D.mep,'nanmean');  % get mean for all subjects/Week/state/task combinations in the entire dataset
[Fs]=pivottable([D.SubjN,D.week,D.state, D.task],[],D.mep,'nanstd');    % get std for all subjects/Week/state/task combinations in the entire dataset

S = [];
for i=1:size(R,1) 
    Di = getrow (D,D.SubjN==R(i,1) & D.week==R(i,2) & D.state==R(i,3) & D.task==R(i,4));
    idx = (Di.mep >F(i)+2*Fs(i)) | (Di.mep <F(i)-2*Fs(i));
    Di = getrow (Di,~idx);
S= addstruct (S,Di);
end;
