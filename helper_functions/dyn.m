D       = dload('bedside_clinical.txt');
D.SID   = strcat(D.Centre,num2str(D.ID));
h       = {'Dyn_FDIL','Dyn_FDIR'};
S       = [];

% 1. Control - non dominant
%   - get dyn measurement for non-dominant hand
T   = getrow(D,D.isPatient==0);
Si  = [];
for i=1:length(T.SID)
    fName       = h{3-T.handedness(i)};
    T.dyn(i,1)  = T.(fName)(i);
end;
Si.SID  = T.SID;
Si.week = T.week;
Si.dyn  = T.dyn;
Si.handLabel = repmat(1,length(T.SID),1);
S = addstruct(S,Si);

% 2. Control - dominant
%   - get dyn measurement for dominant hand
T   = getrow(D,D.isPatient==0);
Si  = [];
for i=1:length(T.SID)
    fName       = h{T.handedness(i)};
    T.dyn(i,1)  = T.(fName)(i);
end;
Si.SID  = T.SID;
Si.week = T.week;
Si.dyn  = T.dyn;
Si.handLabel = repmat(2,length(T.SID),1);
S = addstruct(S,Si);

% 3. Patient - non-paretic
%   - get dyn measurement for non-paretic hand
T   = getrow(D,D.isPatient==1);
Si  = [];
for i=1:length(T.SID)
    fName       = h{3-T.pareticside(i)};
    T.dyn(i,1)  = T.(fName)(i);
end;
Si.SID  = T.SID;
Si.week = T.week;
Si.dyn  = T.dyn;
Si.handLabel = repmat(3,length(T.SID),1);
S = addstruct(S,Si);

% 4. Patient - paretic
%   - get dyn measurement for paretic hand
T   = getrow(D,D.isPatient==1);
Si  = [];
for i=1:length(T.SID)
    fName       = h{T.pareticside(i)};
    T.dyn(i,1)  = T.(fName)(i);
end;
Si.SID  = T.SID;
Si.week = T.week;
Si.dyn  = T.dyn;
Si.handLabel = repmat(4,length(T.SID),1);
S = addstruct(S,Si);


% 5. Plots
figure;
subplot(121);
myboxplot(S.handLabel,S.dyn);
set(gca,'XTickLabel',{'non-dom','dom','non-paretic','paretic'});
ylabel('Force');

subplot(122);
lineplot(S.week,S.dyn,'split',S.handLabel,'plotfcn','nanmean','style_thickline',...
         'leg',{'non-dom','dom','non-paretic','paretic'});
xlabel('week');
ylabel('Force');

keyboard;
