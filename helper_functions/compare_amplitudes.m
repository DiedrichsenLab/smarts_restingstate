function varargout = compare_amplitudes(what,varargin)

%compares IHI state 1 amplitudes with aCST and rCST
%last changed 03/28/2015

D = load ('alldat.mat');
D = getrow (D, D.isGood==1);

% threshold to be exceeded for subject to be considered
count = 9;

[~,R]=pivottable([D.SubjN,D.week],[],D.mep,'length(x)');   % Get all subjects/Week combinations in the entire data set
[F1,R1]=pivottable([D.SubjN,D.week],[],D.mep,'length(x)','subset',D.state==1); % Get all subject/week counts for state 1 & 2 separately
[F2,R2]=pivottable([D.SubjN,D.week],[],D.mep,'length(x)','subset',D.state==2);
S = [];
for i=1:size(R,1)   
    Di = getrow(D,D.SubjN==R(i,1) & D.week==R(i,2));       % get data for subject/week

    % get counts for that subject for states 1 & 2 separately

idx1 = find(R(i,1)==R1(:,1) & R(i,2)==R1(:,2));
idx2 = find(R(i,1)==R1(:,1) & R(i,2)==R1(:,2));
    
    % if either data count does not exist, exclude that subject/week
    try

c1 = F1(idx1);
c2 = F2(idx2);

        if ( c1>count && c2>count )                   % if both states have adequate measurements
        S = addstruct(S,Di);                                % add subject data to new variable S        
        end;        
        catch           % do nothing if there is no measurement
        end;
end;

%for patients
T1 = tapply (S, {'SubjN', 'week', 'task'}, {'mep', 'nanmean', 'subset', S.control==0 & S.state==1 & S.task==2 & S.tms_schedule==0.2000});
T2 = tapply (S, {'SubjN', 'week', 'task'}, {'mep', 'nanmean', 'subset', S.control==0 & S.state==1 & S.task==2 & S.tms_schedule==0.5000});
T3 = tapply (S, {'SubjN', 'week', 'task'}, {'mep', 'nanmean', 'subset', S.control==0 & S.state==1 & S.task==2 & S.tms_schedule==0.8000});
T4 = tapply (S, {'SubjN', 'week', 'task'}, {'mep', 'nanmean', 'subset', S.control==0 & S.state==1 & S.task==2 & S.tms_schedule==0.9500});

T6 = tapply (S, {'SubjN', 'week','task'}, {'rCST', 'nanmean', 'subset', S.control==0 & S.task==2});
T7 = tapply (S, {'SubjN', 'week','task'}, {'aCST', 'nanmean', 'subset', S.control==0 & S.task==2});

% for controls
%T1 = tapply (S, {'SubjN', 'week', 'task'}, {'mep', 'nanmean', 'subset', S.control==1 & S.state==1 & S.task==2 & S.tms_schedule==0.2000});
%T2 = tapply (S, {'SubjN', 'week', 'task'}, {'mep', 'nanmean', 'subset', S.control==1 & S.state==1 & S.task==2 & S.tms_schedule==0.5000});
%T3 = tapply (S, {'SubjN', 'week', 'task'}, {'mep', 'nanmean', 'subset', S.control==1 & S.state==1 & S.task==2 & S.tms_schedule==0.8000});
%T4 = tapply (S, {'SubjN', 'week', 'task'}, {'mep', 'nanmean', 'subset', S.control==1 & S.state==1 & S.task==2 & S.tms_schedule==0.9500});

%T6 = tapply (S, {'SubjN', 'week','task'}, {'rCST', 'nanmean', 'subset', S.control==1 & S.task==2});
%T7 = tapply (S, {'SubjN', 'week','task'}, {'aCST', 'nanmean', 'subset', S.control==1 & S.task==2});

figure; 
barplot (T1.SubjN, [T1.mep, T2.mep, T3.mep, T4.mep, T6.rCST, T7.aCST], 'subset', T1.week==1);

figure;   
barplot (T1.SubjN, [T1.mep, T2.mep, T3.mep, T4.mep, T6.rCST, T7.aCST], 'subset', T1.week==4);

figure;
barplot (T1.SubjN, [T1.mep, T2.mep, T3.mep, T4.mep, T6.rCST, T7.aCST], 'subset', T1.week==12);

figure;
barplot (T1.SubjN, [T1.mep, T2.mep, T3.mep, T4.mep, T6.rCST, T7.aCST], 'subset', T1.week==24);

figure;
barplot (T1.SubjN, [T1.mep, T2.mep, T3.mep, T4.mep, T6.rCST, T7.aCST], 'subset', T1.week==52);

end