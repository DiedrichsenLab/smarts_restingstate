function varargout = get_IHIR(what,varargin)

% get IHIR MEP ratio
% last changed Meret 03/29/16

D = load ('alldat.mat');
D = getrow(D,D.isGood==1 & D.task==1);

% threshold to be exceeded for subject to be considered
count = 9;

[~,R]=pivottable([D.ID,D.week],[],D.mep,'length(x)');   % Get all subjects/Week combinations in the entire data set


[F1,R1]=pivottable([D.ID,D.week],[],D.mep,'length(x)','subset',D.state==1); % Get all subject/week counts for state 1 & 2 separately

[F2,R2]=pivottable([D.ID,D.week],[],D.mep,'length(x)','subset',D.state==2);

 
S = [];

for i=1:size(R,1)                                       % loop over all subject/week measurements

    Di = getrow(D,D.ID==R(i,1) & D.week==R(i,2));       % get data for subject/week

    

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
 
T1 = tapply (S, {'ID', 'week'}, {'mep', 'nanmean', 'subset', S.control==0 & S.state==1});
T2 = tapply (S, {'ID', 'week'}, {'mep', 'nanmean', 'subset', S.control==0 & S.state==2});

T1.quod = T2.mep ./ T1.mep;

myboxplot(T1.week,T1.quod)

% w = unique(T1.week);
% figure
% for i = 1:length(w)  
% subplot(1,length(w),i);
% myboxplot (T1.SubjN, T1.quod,'subset', T1.week==w(i))
% end

