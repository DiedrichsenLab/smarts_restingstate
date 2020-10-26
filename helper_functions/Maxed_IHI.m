function varargout = Maxed_IHI(what,varargin)

%compares IHI state 1 amplitudes with aCST and rCST
%last changed 04/7/2015

D = load ('alldat.mat');
D = getrow (D, D.isGood==1 & D.task==2);

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


T1 = tapply (S, {'SubjN', 'week', 'task'}, {'mep', 'nanmean', 'subset', S.control==0 & S.state==1 & S.tms_schedule==0.9500});
T1a = tapply (S, {'SubjN', 'week', 'task'}, {'mep', 'nanmean', 'subset', S.control==0 & S.state==1 & S.tms_schedule==0.2000});
T2 = tapply (S, {'SubjN', 'week','task'}, {'aCST', 'nanmean', 'subset', S.control==0});

T3 = tapply (S, {'SubjN', 'week', 'task'}, {'mep', 'nanmean', 'subset', S.control==1 & S.state==1 & S.tms_schedule==0.9500});
T3a = tapply (S, {'SubjN', 'week', 'task'}, {'mep', 'nanmean', 'subset', S.control==1 & S.state==1 & S.tms_schedule==0.2000});
T4 = tapply (S, {'SubjN', 'week','task'}, {'aCST', 'nanmean', 'subset', S.control==1});

T1.rel_amp_h =  T1.mep ./ T2.aCST;  
T1a.rel_amp_h = T1a.mep ./ T2.aCST;
T3.rel_amp_h = T3.mep ./ T4.aCST;
T3a.rel_amp_h = T3a.mep ./ T4.aCST;

%patients
P1 = tapply (S, {'SubjN','week', 'tms_schedule'}, {'mep', 'nanmean', 'subset', S.control==0 & S.state==1 & S.tms_schedule==0.2});
P2 = tapply (S, {'SubjN','week', 'tms_schedule'}, {'mep', 'nanmean', 'subset', S.control==0 & S.state==2 & S.tms_schedule==0.2});

P1.IHIP_quod_02 = P2.mep ./ P1.mep; % higher values less inhibition

P3 = tapply (S, {'SubjN','week', 'tms_schedule'}, {'mep', 'nanmean', 'subset', S.control==0 & S.state==1 & S.tms_schedule==0.95});
P4 = tapply (S, {'SubjN','week', 'tms_schedule'}, {'mep', 'nanmean', 'subset', S.control==0 & S.state==2 & S.tms_schedule==0.95});

P1.IHIP_quod_095 = P4.mep ./ P3.mep; % higher values less inhibition

P1.IHIP_evolve = P1.IHIP_quod_095 ./ P1.IHIP_quod_02 ; % higher values the more release of Inhibition

%controls
C1 = tapply (S, {'SubjN','week', 'tms_schedule'}, {'mep', 'nanmean', 'subset', S.control==1 & S.state==1 & S.tms_schedule==0.2});
C2 = tapply (S, {'SubjN','week', 'tms_schedule'}, {'mep', 'nanmean', 'subset', S.control==1 & S.state==2 & S.tms_schedule==0.2});

C1.IHIP_quod_02 = C2.mep ./ C1.mep; % higher values less inhibition

C3 = tapply (S, {'SubjN','week', 'tms_schedule'}, {'mep', 'nanmean', 'subset', S.control==1 & S.state==1 & S.tms_schedule==0.95});
C4 = tapply (S, {'SubjN','week', 'tms_schedule'}, {'mep', 'nanmean', 'subset', S.control==1 & S.state==2 & S.tms_schedule==0.95});

C1.IHIP_quod_095 = C4.mep ./ C3.mep; % higher values less inhibition

C1.IHIP_evolve = C1.IHIP_quod_095 ./ C1.IHIP_quod_02 ; % higher values the more release of Inhibition

figure;
scatterplot(T1.rel_amp_h, P1.IHIP_evolve, 'regression','linear','printcorr', 'split', P1.week, 'leg',{'w1','w4','w12','w24','w52'});
title ('patients');
xlabel ('relative Amplitude state 1 IHIA 0.9/aCST');
ylabel ('IHIA 0.9/0.2 higher values more release of Inhibition');

figure;
scatterplot(T1a.rel_amp_h, P1.IHIP_evolve, 'regression','linear','printcorr', 'split', P1.week, 'leg',{'w1','w4','w12','w24','w52'});
title ('patients'); 
xlabel ('relative Amplitude state 1 IHIA 0.2/aCST');
ylabel ('IHIA 0.9/0.2 higher values more release of Inhibition');

figure;
scatterplot(T3.rel_amp_h, C1.IHIP_evolve, 'regression','linear','printcorr', 'split', C1.week, 'leg',{'w1','w4','w12','w24','w52'});
title ('controls'); 
xlabel ('relative Amplitude state 1 IHIA 0.9/aCST');
ylabel ('IHIA 0.9/0.2 higher values more release of Inhibition');

figure;
scatterplot(T3a.rel_amp_h, C1.IHIP_evolve, 'regression','linear','printcorr', 'split', C1.week, 'leg',{'w1','w4','w12','w24','w52'});
title ('controls'); 
xlabel ('relative Amplitude state 1 IHIA 0.2/aCST');
ylabel ('IHIA 0.9/0.2 higher values more release of Inhibition');




