% get a ratioIHI variable

D = load ('alldat.mat');
D = getrow (D, D.isGood==1 & D.task==2);

for D.mep_outlier == 0
    
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

 S1 = tapply (S, {'subj_name','week', 'tms_schedule'}, {'mep', 'nanmean', 'subset', S.state==1 & S.tms_schedule==0.2});
 S2 = tapply (S, {'subj_name','week', 'tms_schedule'}, {'mep', 'nanmean', 'subset', S.state==2 & S.tms_schedule==0.2});
 S3 = tapply (S, {'subj_name','week', 'tms_schedule'}, {'mep', 'nanmean', 'subset', S.state==1 & S.tms_schedule==0.5});
 S4 = tapply (S, {'subj_name','week', 'tms_schedule'}, {'mep', 'nanmean', 'subset', S.state==2 & S.tms_schedule==0.5});
 S5 = tapply (S, {'subj_name','week', 'tms_schedule'}, {'mep', 'nanmean', 'subset', S.state==1 & S.tms_schedule==0.8});
 S6 = tapply (S, {'subj_name','week', 'tms_schedule'}, {'mep', 'nanmean', 'subset', S.state==2 & S.tms_schedule==0.8});
 S7 = tapply (S, {'subj_name','week', 'tms_schedule'}, {'mep', 'nanmean', 'subset', S.state==1 & S.tms_schedule==0.95});
 S8 = tapply (S, {'subj_name','week', 'tms_schedule'}, {'mep', 'nanmean', 'subset', S.state==2 & S.tms_schedule==0.95});
 
 S1.IHI_02 = S2.mep ./ S1.mep;
 S1.IHI_05 = S4.mep ./ S3.mep;
 S1.IHI_08 = S6.mep ./ S5.mep;
 S1.IHI_095 = S8.mep ./ S7.mep;
 
 
 %save as
 
 