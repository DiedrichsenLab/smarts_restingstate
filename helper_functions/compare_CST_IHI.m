function varargout = compare_CST_IHI(what,varargin)

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
for i=1:size(R,1)                                       % loop over all subject/week measurements
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


T1 = tapply (S, {'subj_name', 'week', 'task'}, {'mep', 'nanmean', 'subset', S.control==0 & S.state==1 & S.task==2 & S.tms_schedule==0.9500});
T1a = tapply (S, {'subj_name', 'week', 'task'}, {'mep', 'nanmean', 'subset', S.control==0 & S.state==1 & S.task==2 & S.tms_schedule==0.2000});
T2 = tapply (S, {'subj_name', 'week','task'}, {'aCST', 'nanmean', 'subset', S.control==0 & S.task==2});

T3 = tapply (S, {'subj_name', 'week', 'task'}, {'mep', 'nanmean', 'subset', S.control==1 & S.state==1 & S.task==2 & S.tms_schedule==0.9500});
T3a = tapply (S, {'subj_name', 'week', 'task'}, {'mep', 'nanmean', 'subset', S.control==1 & S.state==1 & S.task==2 & S.tms_schedule==0.2000});
T4 = tapply (S, {'subj_name', 'week','task'}, {'aCST', 'nanmean', 'subset', S.control==1 & S.task==2});

T1.rel_amp_h = T1.mep ./ T2.aCST;    % for last TMS_IHIP_timing 0.95
T1a.rel_amp_h = T1a.mep ./ T2.aCST;  % for last TMS_IHIP_timing 0.2
T3.rel_amp_h = T3.mep ./ T4.aCST ;   % for last TMS_IHIP_timing 0.95
T3a.rel_amp_h =T3a.mep ./ T4.aCST ;  % for last TMS_IHIP_timing 0.2
    
figure;
scatterplot(T1.week,T1.rel_amp_h);
title ('patients')
xlabel ('Week');
ylabel ('IHIA TS  0.95/ aCST (mV)');

figure;
scatterplot(T3.week,T3.rel_amp_h);
title ('controls');
xlabel ('Week');
ylabel ('IHIA TS  0.95/ aCST (mV)');

figure;
scatterplot(T1a.week,T1a.rel_amp_h);
title ('patients');
xlabel ('Week');
ylabel ('IHIA TS  0.2/ aCST (mV)');

figure;
scatterplot(T3a.week,T3a.rel_amp_h);
title ('controls');
xlabel ('Week');
ylabel ('IHIA TS  0.2/ aCST (mV)');



