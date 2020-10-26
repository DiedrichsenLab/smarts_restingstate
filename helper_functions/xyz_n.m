% IHI_clinics

D = load ('alldat.mat');
D = getrow (D, D.isGood==1 & D.task==2);%

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

% 1. Week/Schedule plot for IHI
S = tapply (S, {'SubjN','week', 'tms_schedule', 'control','state'}, {'mep', 'nanmean'});
S1 = getrow(S,S.state==1);
S2 = getrow(S,S.state==2);

S1.IHI = S2.mep./S1.mep;
S = rmfield(S1,{'mep','state'});

tm      = unique(S.tms_schedule);
week    = unique(S.week);


for i=1:length(week)
    figure;
    for j=1:length(tm)
        
        Si = getrow(S,S.week==week(i) & S.tms_schedule==tm(j));
        
        %Plotting each individual week/tms schedule
        subplot(1,4,j);
        myboxplot(Si.control,Si.IHI);
        set(gca,'XTickLabel',{'patients','controls'});
        ylabel('IHI');
        ylim([0 2.5]);
        
        
        
        % statistics
        fprintf('Week %d, Timing %2.2f\n',week(i),tm(j));
        fprintf('------------------------\n');
        ttest(Si.IHI(Si.control==0),Si.IHI(Si.control==1),2,'independent');
        fprintf('\n\n');
    end;
end;
      

% Check the tms schedule
pivottable([S.SubjN, S.week],S.tms_schedule,S.IHI,...
           'nanmean','subset',ismember(S.tms_schedule,[0.2 0.95]) & S.control==0);
       
pivottable([S.SubjN, S.week],S.tms_schedule,S.IHI,...
           'nanmean','subset',ismember(S.tms_schedule,[0.2 0.95]) & S.control==1);

% 2. Delta IHI plot for Week


dcount=0;

[~,Rd]=pivottable([S.SubjN, S.week],[],S.IHI,'length(x)');   % Get all subjects/Week combinations in the entire data set
[F1d,R1d]=pivottable([S.SubjN, S.week],[],S.IHI,'length(x)','subset',S.tms_schedule==0.2); % Get all subject/week counts for state 1 & 2 separately
[F2d,R2d]=pivottable([S.SubjN, S.week],[],S.IHI,'length(x)','subset',S.tms_schedule==0.95);
Sd = [];
for id=1:size(Rd,1)                                      % loop over all subject/week measurements
Did = getrow(S,S.SubjN==Rd(id,1) & S.week==Rd(id,2));    % get data for subject/week
% get counts for that subject for states 1 & 2 separately
idx1d = find(Rd(id,1)==R1d(:,1) & Rd(id,2)==R1d(:,2));
idx2d = find(Rd(id,1)==R2d(:,1) & Rd(id,2)==R2d(:,2));
% if either data count does not exist, exclude that subject/week
try
cd1 = F1d(idx1d);
cd2 = F2d(idx2d);
if ( cd1>dcount && cd2>dcount )                   % if both states have adequate measurements
Sd = addstruct(Sd,Did);                                % add subject data to new variable S
end;
catch           % do nothing if there is no measurement
end;
end

T1 = getrow(Sd,Sd.tms_schedule==0.2);
T2 = getrow(Sd,Sd.tms_schedule==0.95);

T1.delta_IHI = T2.IHI./T1.IHI;

T = rmfield(T1,{'tms_schedule','IHI'});

week    = unique(T.week);

figure;
for i=1:length(week)
    Ti = getrow(T,T.week==week(i));

    % Plotting each individual week/tms schedule
    subplot(1,5,i);
    myboxplot(Ti.control,Ti.delta_IHI);
    set(gca,'XTickLabel',{'patients','controls'});
    ylabel('IHI');
    ylim([0 2.5]);
    title(sprintf('Week %d\n',week(i)));

    % statistics
    fprintf('Week %d\n',week(i));
    fprintf('--------\n');
    ttest(Ti.delta_IHI(Ti.control==0),Ti.delta_IHI(Ti.control==1),2,'independent');
    fprintf('\n\n');
end;

% 3. Statistics
%   - average over weeks
%   - anova for controls
%   - anova for patients
X = tapply(S,{'SubjN','tms_schedule','control'},{'IHI','nanmean'});

F=pivottable(X.SubjN,S.tms_schedule,X.IHI,'nanmean','subset',X.control==1);
anova1(F);

F=pivottable(X.SubjN,S.tms_schedule,X.IHI,'nanmean','subset',X.control==0);
anova1(F);

%  S1 = tapply (S, {'subj_name','week', 'tms_schedule', 'control'}, {'mep', 'nanmean', 'subset', S.state==1 & S.tms_schedule==0.2});
%  S2 = tapply (S, {'subj_name','week', 'tms_schedule', 'control'}, {'mep', 'nanmean', 'subset', S.state==2 & S.tms_schedule==0.2});
%  S3 = tapply (S, {'subj_name','week', 'tms_schedule', 'control'}, {'mep', 'nanmean', 'subset', S.state==1 & S.tms_schedule==0.5});
%  S4 = tapply (S, {'subj_name','week', 'tms_schedule', 'control'}, {'mep', 'nanmean', 'subset', S.state==2 & S.tms_schedule==0.5});
%  S5 = tapply (S, {'subj_name','week', 'tms_schedule', 'control'}, {'mep', 'nanmean', 'subset', S.state==1 & S.tms_schedule==0.8});
%  S6 = tapply (S, {'subj_name','week', 'tms_schedule', 'control'}, {'mep', 'nanmean', 'subset', S.state==2 & S.tms_schedule==0.8});
%  S7 = tapply (S, {'subj_name','week', 'tms_schedule', 'control'}, {'mep', 'nanmean', 'subset', S.state==1 & S.tms_schedule==0.95});
%  S8 = tapply (S, {'subj_name','week', 'tms_schedule', 'control'}, {'mep', 'nanmean', 'subset', S.state==2 & S.tms_schedule==0.95});
%  
%  S1.IHI_02 = S2.mep ./ S1.mep;
%  S1.IHI_05 = S4.mep ./ S3.mep;
%  S1.IHI_08 = S6.mep ./ S5.mep;
%  S1.IHI_095 = S8.mep ./ S7.mep;
% 
%  
% % make state 1 comparing graphs
% figure;
% barplot (S.control, S.mep, 'subset', S.state==1 & S.week==1 & S.tms_schedule==0.2);
% ylabel ('TS alone MEP mV')
% title ('Week 1 0.2')
% figure;
% barplot (S.control, S.mep, 'subset', S.state==1 & S.week==1 & S.tms_schedule==0.5);
% ylabel ('TS alone MEP mV')
% title ('Week 1 0.5')
% figure;
% barplot (S.control, S.mep, 'subset', S.state==1 & S.week==1 & S.tms_schedule==0.8);
% ylabel ('TS alone MEP mV')
% title ('Week 1 0.8')
% figure;
% barplot (S.control, S.mep, 'subset', S.state==1 & S.week==1 & S.tms_schedule==0.95);
% ylabel ('TS alone MEP mV')
% title ('Week 1 0.9')
% 
% % IHI data plotting
% weeks = [1 4 12 24 52];
% for w=1:length(weeks);
%     
% figure;
% subplot(141);
% myboxplot(S1.control,S1.IHI_02, 'subset', S1.week==weeks(w));
% set(gca,'XTickLabel',{'patients','controls'});
% ylabel('IHI');
% ylim([0 2.5]);
% 
% subplot(142);
% myboxplot(S1.control,S1.IHI_05, 'subset', S1.week==weeks(w));
% set(gca,'XTickLabel',{'patients','controls'});
% ylabel('IHI');
% ylim([0 2.5]);
% 
% subplot(143);
% myboxplot(S1.control,S1.IHI_08, 'subset', S1.week==weeks(w));
% set(gca,'XTickLabel',{'patients','controls'});
% ylabel('IHI');
% ylim([0 2.5]);
% 
% subplot(144);
% myboxplot(S1.control,S1.IHI_095, 'subset', S1.week==weeks(w));
% set(gca,'XTickLabel',{'patients','controls'});
% ylabel('IHI');
% ylim([0 2.5]);
% 
% end;
% 
% % ttest
% 
% var = {'IHI_02','IHI_05','IHI_08','IHI_095'};
% weeks = [1 4 12 24 52];
% 
% for w=1:length(weeks)
%     for v=1:length(var)
%         fprintf('Week: %d, TP: %s\n-----------------------\n',weeks(w),var{v});
%         F = pivottable(S1.subj_name,S1.control,S1.(var{v}),'nanmean','subset',S1.week==weeks(w));
%         ttest(F(:,1),F(:,2),2,'independent');
%         fprintf('\n\n');
%     end;
% end;
% 
% 
% 
% S1a = tapply (S, {'subj_name','week', 'tms_schedule', 'control', 'FM_hand'}, {'mep', 'nanmean', 'subset', D.state==1 & D.tms_schedule==0.95 & S.control==0 & S.week==1});
%  S2a = tapply (S, {'subj_name','week', 'tms_schedule', 'control', 'FM_hand'}, {'mep', 'nanmean', 'subset', S.state==2 & S.tms_schedule==0.95 & S.control==0 & S.week==1});
%  S3a = tapply (S, {'subj_name','week', 'tms_schedule', 'control', 'FM_hand'}, {'mep', 'nanmean', 'subset', S.state==1 & S.tms_schedule==0.95 & S.control==0 & S.week==4});
%  S4a = tapply (S, {'subj_name','week', 'tms_schedule', 'control', 'FM_hand'}, {'mep', 'nanmean', 'subset', S.state==2 & S.tms_schedule==0.95 & S.control==0 & S.week==4});
%  S5a = tapply (S, {'subj_name','week', 'tms_schedule', 'control', 'FM_hand'}, {'mep', 'nanmean', 'subset', S.state==1 & S.tms_schedule==0.95 & S.control==0 & S.week==12});
%  S6a = tapply (S, {'subj_name','week', 'tms_schedule', 'control', 'FM_hand'}, {'mep', 'nanmean', 'subset', S.state==2 & S.tms_schedule==0.95 & S.control==0 & S.week==12});
%  S7a = tapply (S, {'subj_name','week', 'tms_schedule', 'control', 'FM_hand'}, {'mep', 'nanmean', 'subset', S.state==1 & S.tms_schedule==0.95 & S.control==0 & S.week==24});
%  S8a = tapply (S, {'subj_name','week', 'tms_schedule', 'control', 'FM_hand'}, {'mep', 'nanmean', 'subset', S.state==2 & S.tms_schedule==0.95 & S.control==0 & S.week==24});
%  S9a = tapply (S, {'subj_name','week', 'tms_schedule', 'control', 'FM_hand'}, {'mep', 'nanmean', 'subset', S.state==1 & S.tms_schedule==0.95 & S.control==0 & S.week==52});
%  S10a = tapply (S, {'subj_name','week', 'tms_schedule', 'control', 'FM_hand'}, {'mep', 'nanmean', 'subset', S.state==2 & S.tms_schedule==0.95 & S.control==0 & S.week==52});
% 
%  S1a.IHI_week1 = S2a.mep ./ S1a.mep;
% S4a.IHI_week2 = S4a.mep ./ S3a.mep;
% S6a.IHI_week3 = S6a.mep ./ S5a.mep;
% S8a.IHI_week4 = S8a.mep ./ S7a.mep;
%  
% scatterplot( S1a.IHI_week1, S1a.FM_hand, 'regression', 'linear', 'split',S1a.subj_name,'subset', S1a.control==0);
% figure;
% scatterplot( S4a.IHI_week2, S4a.FM_hand, 'regression', 'linear', 'split',S4a.subj_name,'subset', S4a.control==0);
% figure;
% scatterplot( S6a.IHI_week3, S6a.FM_hand, 'regression', 'linear', 'split',S6a.subj_name,'subset', S6a.control==0);
% figure;
% scatterplot( S8a.IHI_week4, S8a.FM_hand, 'regression', 'linear', 'split',S8a.subj_name,'subset', S8a.control==0);