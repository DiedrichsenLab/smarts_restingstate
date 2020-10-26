%T-testing for IHIR

D = load ('alldat.mat');
D = getrow (D, D.isGood==1 & D.task==1);

% S=getSubjWeek(D,1); add if you only want to have subjects with >1 time
% point

count=9;

[~,R]=pivottable([D.SubjN, D.week, D.tms_schedule],[],D.mep,'length(x)');   % Get all subjects/Week combinations in the entire data set
[F1,R1]=pivottable([D.SubjN,D.week, D.tms_schedule],[],D.mep,'length(x)','subset',D.state==1); % Get all subject/week counts for state 1 & 2 separately
[F2,R2]=pivottable([D.SubjN,D.week, D.tms_schedule],[],D.mep,'length(x)','subset',D.state==2);
S = [];
for i=1:size(R,1)                                      % loop over all subject/week measurements
Di = getrow(D,D.SubjN==R(i,1) & D.week==R(i,2) & D.tms_schedule==R(i,3));       % get data for subject/week
% get counts for that subject for states 1 & 2 separately
idx1 = find(R(i,1)==R1(:,1) & R(i,2)==R1(:,2) & R(i,3)==R1(:,3));
idx2 = find(R(i,1)==R2(:,1) & R(i,2)==R2(:,2) & R(i,3)==R2(:,3));
% if either data count does not exist, exclude that subject/week
try
c1 = F1(idx1);
c2 = F2(idx2);
if ( c1>count && c2>count )                   % if both states have adequate measurements
S = addstruct(S,Di);                                % add subject data to new variable S
end;
catch           % do nothing if there is no measurement
end;
end
D = getrow(D,D.task==1);
tests = [1 1; 1 4; 1 12; 1 24; 1 52; 2 1; 2 4; 2 12; 2 24; 2 52];

% MEP (mV)
% for i = 1:size(tests,1)
% ttestDirect(S.mep,[S.control S.SubjN],2,'independent','subset',S.state==tests(i,1) & S.week==tests(i,2));
% end;

% (CS+TS)/T
% T1 = tapply (S, {'SubjN', 'week', 'control'}, {'mep', 'nanmean', 'subset', S.state==1});
% T2 = tapply (S, {'SubjN', 'week', 'control'}, {'mep', 'nanmean', 'subset', S.state==2});

% T1.quod = T2.mep ./ T1.mep;



% for i = 1:size(tests,1)
% ttestDirect(T1.quod,[T1.control T1.SubjN],2,'independent','subset', T1.week==tests(i,2));
% end;

%careful with intensities because CS_intensity = state 1 and state 2
%only five last rows
for i = 1:size(tests,1)
ttestDirect(S.CS_intensity,[S.control S.SubjN],2,'independent','subset',S.state==tests(i,1) & S.week==tests(i,2));
end;

% for i = 1:size(tests,1)
% ttestDirect(S.TS_intensity,[S.control S.SubjN],2,'independent','subset',S.state==tests(i,1) & S.week==tests(i,2));
% end;