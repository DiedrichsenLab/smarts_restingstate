function varargout=rs_imana_new(what,varargin)
% MB - analysis script for the SMARTS resting state data
% Open questions:
% - r or fisher z
% - dominant = nonlesion controls???


% rootDir             = '/Volumes/External/data/smarts';
% rootDir             = '/Users/naveed/Documents/data/smarts';
% rootDir               = '/Users/meret/Dropbox/N&M/rsfMRI_results';
% rootDir              = '/Users/naveed/Dropbox/Personal/N&M/rsfMRI_results';
rootDir             = '/Volumes/Naveed/data/smarts';
behDir              = [rootDir '/bedside/analysis'];
rsDir               = [rootDir '/fmri/restingstate_imaging'];
roiDir              = [rootDir '/fmri/RegionOfInterest'];
ppDir               = [rsDir '/preprocess/'];

% rootDir             = ['~/Dropbox/Personal/N&M/rsfMRI_results/new'];
rDir                = [rootDir, '/R'];

% color_scheme
colours             = num2cell(colormap(parula(20)),2)';        close;
sty_2grp_ac         = colours([8,12]);   % e.g. group 1 vs group 2
sty_multiple        = colours([10;14;16]);

regSide             = [ones(1,8), ones(1,8)+1]';
regType             = [[1:8]'; [1:8]'];

regLabel = {'S1-S1','S1-M1','S1-PMd','S1-PMv','S1-SMA',...
            'M1-S1','M1-M1','M1-PMd','M1-PMv','M1-SMA',...
            'PMd-S1','PMd-M1','PMd-PMd','PMd-PMv','PMd-SMA',...
            'PMv-S1','PMv-M1','PMv-PMd','PMv-PMv','PMv-SMA',...
            'PMd-S1','PMd-M1','PMd-PMd','PMd-PMv','PMd-SMA'};

switch(what)
    case 'unzip'
        cd(ppDir);
        d = dir('*.gz');
        for i=1:length(d)
            fprintf('%d. %s\n',i,d(i).name);
            try
                gunzip(d(i).name);
            catch
                disp('unable to process');
            end;
        end;
    case 'preprocess_comparison'
        type = {'pre','preproc'};
        roi  = [2 10];  % left and right M1
        
        cd(ppDir);
        d = dir('*.nii');

        S = [];
        for i=1:length(d)
            fprintf('%s\n',d(i).name);
            s   = textscan(d(i).name,'%s%s%s%s%s%s.nii','delimiter','_');
            w   = strtok(s{end},'W.nii');
            sn  = strcat(s{4},'_',s{5});
            
            R   = load(fullfile(roiDir,sprintf('%s_regions_all.mat',sn{1}))); % load subject region definition
            v   = spm_vol(d(i).name);     % load volumes from rs data
            ts  = region_getdata(v,R.R(roi));   % get time series for M1-left and M1-right

            % calculate corr for time series
            ts1 = mean(ts{1},2);
            ts2 = mean(ts{2},2);
            r   = corr(ts1,ts2);
            
            Si.subj_name    = sn;
            Si.week         = str2double(w{1});
            Si.type         = find(strcmp(s{1},type));            
            Si.r            = r;
            Si.rfisherz     = fisherz(r);
            S               = addstruct(S,Si);
        end;
        varargout = {S};
        save(fullfile(ppDir,'rs_preprocess.mat'),'-struct','S');
        
    case 'PP_rawTimeSeries'             % Extraction of fisher-z correlations from different regions of interest from the freesurfer atlas
        pre = {'preproc_Bold_Rest'};
        lesionType = {'none','cortical','subcortical','mixed'}; %1,2,3,4
        roi  = combnk(1:16,2); % all pairs
        
        cd(ppDir);
        d = dir('*.nii');
        
        D   = dload(fullfile(ppDir,'patient_list.txt'));
        B   = load(fullfile(behDir,'ens_alldat_mvc.mat'));
        SS  = dload(fullfile(rootDir,'bedside','data_clean','subject_list.txt'));
        L   = dload(fullfile(rootDir,'DTI','lesion_analysis','lesion.dat'));
        
        SS.subj_name    = strcat(SS.Centre,'_',num2str(SS.ID));
        L.subj_name     = strcat(L.Centre,'_',num2str(L.ID));
        
        S = [];
        for i=1:length(D.SN)
            fname   = strcat(pre{1},'_',D.Centre{i},'_',num2str(D.ID(i)),'_',D.Week{i},'.nii');
            sn      = strcat(D.Centre{i},'_',num2str(D.ID(i)));
            fprintf('%d.\t%s\n',D.SN(i),fname);
            
            regFile = fullfile(roiDir,sprintf('%s_regions_all.mat',sn));
            if ~exist(regFile)
                fprintf('No region definition for this subject...\n');
            else
                R   = load(fullfile(roiDir,sprintf('%s_regions_all.mat',sn))); % load subject region definition
                v   = spm_vol(d(i).name);       % load volumes from rs data
                ts  = region_getdata(v,R.R);    % get time series for all ROIs
                Bi  = getrow(B,strcmp(B.subj_name,sn) & B.week==D.week(i));
                SSi = getrow(SS,strcmp(SS.subj_name,sn));   % get gray matter damage from subject_list file
                Li  = getrow(L,strcmp(L.subj_name,sn));     % get cst damage from DTI
                
                clear Si
                Si.subj_name    = {sn};
                Si.week         = D.week(i);
                Si.lesionType   = find(strcmp(D.LesionLocation(i),lesionType));
                if isempty(find(strcmp(D.LesionLocation(i),lesionType)))
                    Si.lesionType = nan;
                end;
                Si.lesionSide   = D.lesionside(i);                          %lesionside = 1 = left??? 
                Si.control      = ~isempty(strfind(sn,'P'));
                if ~isempty(Bi.subj_name)
                    switch(Si.control)
                        case 0
                            Bj = getrow(Bi,Bi.handLabel==4);    % paretic hand
                        case 1  
                            Bj = getrow(Bi,Bi.handLabel==1);    % non-dom hand
                    end;
                    Si.ensOverall   = -Bj.ensOverall;
                    Si.mmOverall    = -Bj.mmOverall;
                    Si.mvcNorm      = Bj.mvcNorm;
                    Si.FM           = Bj.FM;
                    Si.lesion_PrCG  = SSi.lesion_PrCG(1);
                    Si.lesion_PrCG_W= SSi.lesion_PrCG_W(1);
                else
                    Si.ensOverall   = nan;
                    Si.mmOverall    = nan;
                    Si.mvcNorm      = nan;
                    Si.FM           = nan;
                    Si.lesion_PrCG  = nan;
                    Si.lesion_PrCG_W= nan;
                end;
                if (~isempty(Li.subj_name))
                    Si.lesion_SubCW4_noCST          = Li.lesion_SubCW4_noCST;
                    Si.lesion_SubCW3_noCST          = Li.lesion_SubCW3_noCST;
                    Si.lesion_total_volume_noCST    = Li.lesion_total_volume_noCST;
                    Si.lesion_total_volume          = Li.lesion_total_volume;
                    Si.lesion_total_fa_noCST        = Li.lesion_total_fa_noCST;
                    Si.lesion_total_fa              = Li.lesion_total_fa;
                    Si.lesion_CST_total             = Li.lesion_CST_total;
                    Si.CST_Hemi                     = Li.CST_Hemi;
                    Si.lesion_SubCW7_noCST          = Li.lesion_SubCW7_noCST;
                    Si.lesion_SubCW7                = Li.lesion_SubCW7;
                    Si.lesion_SubCW7_CST            = Li.lesion_SubCW7_CST;
                    Si.lesion_CST_SP2               = Li.lesion_CST_SP2;
                else
                    Si.lesion_SubCW4_noCST          = nan;
                    Si.lesion_SubCW3_noCST          = nan;
                    Si.lesion_total_volume_noCST    = nan;
                    Si.lesion_total_volume          = nan;
                    Si.lesion_total_fa_noCST        = nan;
                    Si.lesion_total_fa              = nan;
                    Si.lesion_CST_total             = nan;
                    Si.CST_Hemi                     = nan;
                    Si.lesion_SubCW7_noCST          = nan;
                    Si.lesion_SubCW7                = nan;
                    Si.lesion_SubCW7_CST            = nan;
                    Si.lesion_CST_SP2               = nan;
                end;

                for r=1:size(roi,1)
                    rpair = roi(r,:);
                    ts1 = nanmean(ts{rpair(1)},2);
                    ts2 = nanmean(ts{rpair(2)},2);
                    rp  = corr(ts1,ts2);
                    
                    Si.rType    = r;
                    Si.rpair    = rpair;
                    Si.BOLD     = [nanmean(ts1) nanmean(ts2)];
                    Si.r        = rp;
                    Si.fz_r     = fisherz(rp);     % 1. save fisher z transform r values                    
                    % 2. add location for region side (left or right hem)
                    %   - 1 is left hem, 2 is right
                    Si.hem      = regSide(rpair)';     
                    % 3. add a variable call hl which is
                    %   - 0 if interhem, 1 if intra (left hem), 2 if intra (right hem)
                    Si.hl       = mean(Si.hem,2);
                    Si.hl(Si.hl==1.5) = 0;  
                    % 4. add a variable called regType which is
                    %   - 1 for S1, 2 for M1, 3 for PMd etc ...
                    Si.regType  = regType(Si.rpair)';                    
                    
                    S           = addstruct(S,Si);
                end;
            end;
        end;
        varargout = {S};
        save(fullfile(ppDir,'rs_preprocess.mat'),'-struct','S');            
    case 'PP_rawTimeSeries_splitHalf'   % Extraction of fisher-z correlations from different regions of interest from the freesurfer atlas
        pre = {'preproc_Bold_Rest'};
        lesionType = {'none','cortical','subcortical','mixed'}; %1,2,3,4
        roi  = combnk(1:16,2); % all pairs
        
        cd(ppDir);
        d = dir('*.nii');
        
        D   = dload(fullfile(ppDir,'patient_list.txt'));
        B   = load(fullfile(behDir,'ens_alldat_mvc.mat'));
        SS  = dload(fullfile(rootDir,'bedside','data_clean','subject_list.txt'));
        L   = dload(fullfile(rootDir,'DTI','lesion_analysis','lesion.dat'));
        
        SS.subj_name    = strcat(SS.Centre,'_',num2str(SS.ID));
        L.subj_name     = strcat(L.Centre,'_',num2str(L.ID));
        
        S = [];
        for i=1:length(D.SN)
            fname   = strcat(pre{1},'_',D.Centre{i},'_',num2str(D.ID(i)),'_',D.Week{i},'.nii');
            sn      = strcat(D.Centre{i},'_',num2str(D.ID(i)));
            fprintf('%d.\t%s\n',D.SN(i),fname);
            
            regFile = fullfile(roiDir,sprintf('%s_regions_all.mat',sn));
            if ~exist(regFile)
                fprintf('No region definition for this subject...\n');
            else
                R   = load(fullfile(roiDir,sprintf('%s_regions_all.mat',sn))); % load subject region definition
                v   = spm_vol(d(i).name);       % load volumes from rs data
                nVol= length(v);                  % total number of volumes
                ts  = region_getdata(v,R.R);    % get time series for all ROIs
                Bi  = getrow(B,strcmp(B.subj_name,sn) & B.week==D.week(i));
                SSi = getrow(SS,strcmp(SS.subj_name,sn));   % get gray matter damage from subject_list file
                Li  = getrow(L,strcmp(L.subj_name,sn));     % get cst damage from DTI
                
                clear Si
                Si.subj_name    = {sn};
                Si.week         = D.week(i);
                Si.lesionType   = find(strcmp(D.LesionLocation(i),lesionType));
                if isempty(find(strcmp(D.LesionLocation(i),lesionType)))
                    Si.lesionType = nan;
                end;
                Si.lesionSide   = D.lesionside(i);                          %lesionside = 1 = left??? 
                Si.control      = ~isempty(strfind(sn,'P'));
                if ~isempty(Bi.subj_name)
                    switch(Si.control)
                        case 0
                            Bj = getrow(Bi,Bi.handLabel==4);    % paretic hand
                        case 1  
                            Bj = getrow(Bi,Bi.handLabel==1);    % non-dom hand
                    end;
                    Si.ensOverall   = -Bj.ensOverall;
                    Si.mmOverall    = -Bj.mmOverall;
                    Si.mvcNorm      = Bj.mvcNorm;
                    Si.FM           = Bj.FM;
                    Si.lesion_PrCG  = SSi.lesion_PrCG(1);
                    Si.lesion_PrCG_W= SSi.lesion_PrCG_W(1);
                else
                    Si.ensOverall   = nan;
                    Si.mmOverall    = nan;
                    Si.mvcNorm      = nan;
                    Si.FM           = nan;
                    Si.lesion_PrCG  = nan;
                    Si.lesion_PrCG_W= nan;
                end;
                if (~isempty(Li.subj_name))
                    Si.lesion_SubCW4_noCST          = Li.lesion_SubCW4_noCST;
                    Si.lesion_SubCW3_noCST          = Li.lesion_SubCW3_noCST;
                    Si.lesion_total_volume_noCST    = Li.lesion_total_volume_noCST;
                    Si.lesion_total_volume          = Li.lesion_total_volume;
                    Si.lesion_total_fa_noCST        = Li.lesion_total_fa_noCST;
                    Si.lesion_total_fa              = Li.lesion_total_fa;
                    Si.lesion_CST_total             = Li.lesion_CST_total;
                    Si.CST_Hemi                     = Li.CST_Hemi;
                    Si.lesion_SubCW7_noCST          = Li.lesion_SubCW7_noCST;
                    Si.lesion_SubCW7                = Li.lesion_SubCW7;
                    Si.lesion_SubCW7_CST            = Li.lesion_SubCW7_CST;
                    Si.lesion_CST_SP2               = Li.lesion_CST_SP2;
                else
                    Si.lesion_SubCW4_noCST          = nan;
                    Si.lesion_SubCW3_noCST          = nan;
                    Si.lesion_total_volume_noCST    = nan;
                    Si.lesion_total_volume          = nan;
                    Si.lesion_total_fa_noCST        = nan;
                    Si.lesion_total_fa              = nan;
                    Si.lesion_CST_total             = nan;
                    Si.CST_Hemi                     = nan;
                    Si.lesion_SubCW7_noCST          = nan;
                    Si.lesion_SubCW7                = nan;
                    Si.lesion_SubCW7_CST            = nan;
                    Si.lesion_CST_SP2               = nan;
                end;

                for r=1:size(roi,1)
                    rpair = roi(r,:);
                    ts1 = nanmean(ts{rpair(1)},2);
                    ts2 = nanmean(ts{rpair(2)},2);
                    
                    nSplit = {1:floor(nVol/2),floor(nVol/2)+1:nVol};
                    for n=1:length(nSplit)
                        split1      = ts1(nSplit{n});
                        split2      = ts2(nSplit{n});
                        
                        rp          = corr(split1,split2);  % split half correlations                        
                        
                        Si.rType    = r;
                        Si.rpair    = rpair;
                        Si.BOLD     = [nanmean(split1) nanmean(split2)];
                        Si.split    = n;                        
                        Si.r        = rp;

                        Si.fz_r     = fisherz(rp);     % 1. save fisher z transform r values                    
                        % 2. add location for region side (left or right hem)
                        %   - 1 is left hem, 2 is right
                        Si.hem      = regSide(rpair)';     
                        % 3. add a variable call hl which is
                        %   - 0 if interhem, 1 if intra (left hem), 2 if intra (right hem)
                        Si.hl       = mean(Si.hem,2);
                        Si.hl(Si.hl==1.5) = 0;  
                        % 4. add a variable called regType which is
                        %   - 1 for S1, 2 for M1, 3 for PMd etc ...
                        Si.regType  = regType(Si.rpair)';

                        S           = addstruct(S,Si);
                    end;
                end;
            end;
        end;
        varargout = {S};
        save(fullfile(ppDir,'rs_preprocess_splithalf.mat'),'-struct','S');                
    case 'PP_removeInf'                 % some participants have infinity fisherz correlations, remove these time points
        D = varargin{1};
        
        % 0. Figure out which participants have inf in their fisherz scores
        [x,sn,week]     = pivottable(D.subj_name,D.week,D.fz_r,'sum(isinf(x))');
        [i,j]           = ind2sub(size(x),find(x>0));
        
        idx = zeros(length(D.subj_name),1);    
        for k=1:length(i)
            l       = strcmp(D.subj_name,sn(i(k))) & D.week==week(j(k));
            idx(l)  = 1;
            fprintf('%s, week %d removed\n',sn{i(k)},week(j(k)));
        end;
        D = getrow(D,~idx);
        varargout = {D};
        
    case 'NC_getInterHemPattern'           % get pattern of correlations for pairs between hemispheres
        % 0. Load data
        D       = load(fullfile(rootDir,'rs_preprocess.mat')); 
        D       = rs_imana_new('PP_removeInf',D);

        % 0. Obtain required data subset        
        %       - exclude unwanted ROIs
        inclROI = [1:5];        %no parietal or v1 rois         
        idx = ismember(D.regType(:,1),inclROI) & ismember(D.regType(:,2),inclROI);        
        D   = getrow(D,idx);        
        
        %       - get interhemispheric fisherz correlations
        D   = getrow(D,D.hl==0);
        disp(what);
        
        S   = [];
        % Loop over subj
        for sn=unique(D.subj_name)'
            Ds = getrow(D,strcmp(D.subj_name,sn));
            
            % Loop over week
            for w=unique(Ds.week)'
                Dw  = getrow(Ds,Ds.week==w);
                fprintf('%s - W%d\n',Dw.subj_name{1},Dw.week(1));
                
                M   = zeros(length(inclROI),length(inclROI));     % matrix of corr. for all regions
                for i=1:size(Dw.regType,1)
                    M(Dw.regType(i,1),Dw.regType(i,2)) = Dw.fz_r(i);
                end;
                
                % add data to new structure
                Si      = getrow(Dw,1);
                Si      = rmfield(Si,{'rType','rpair','r'});
                Si.M    = nonsym_squareform(M);
                S       = addstruct(S,Si);
            end;
        end;
        save(fullfile(rootDir,'nc_interhem.mat'),'-struct','S'); 
    case 'NC_getIntraLesionedHemPattern'   % get pattern of correlations for pairs within lesioned hemisphere
        % 0. Load data
        D       = load(fullfile(rootDir,'rs_preprocess.mat')); 
        D       = rs_imana_new('PP_removeInf',D);

        % 0. Obtain required data subset        
        %       - exclude unwanted ROIs
        inclROI = [1:5];        %no parietal or v1 rois         
        idx = ismember(D.regType(:,1),inclROI) & ismember(D.regType(:,2),inclROI);        
        D   = getrow(D,idx);        
        
        %       - get intrahemispheric fisherz correlations
        D   = getrow(D,D.hl~=0);
        Dc  = getrow(D,D.control==1);
        Dp  = getrow(D,D.control==0 & D.lesionSide==D.hl);
        D   = addstruct(Dp,Dc);
        disp(what);
        
        S   = [];
        % Loop over subj
        for sn=unique(D.subj_name)'
            Ds = getrow(D,strcmp(D.subj_name,sn));
            
            % Loop over week
            for w=unique(Ds.week)'
                Dw  = getrow(Ds,Ds.week==w);
                fprintf('%s - W%d\n',Dw.subj_name{1},Dw.week(1));
                
                M   = zeros(length(inclROI),length(inclROI));     % matrix of corr. for all regions
                for i=1:size(Dw.regType,1)
                    M(Dw.regType(i,1),Dw.regType(i,2)) = Dw.fz_r(i);
                end;
                M   = M+M';
                
                % add data to new structure
                Si      = getrow(Dw,1);
                Si      = rmfield(Si,{'rType','rpair','r'});
                Si.M    = squareform(M);
                S       = addstruct(S,Si);
            end;
        end;
        save(fullfile(rootDir,'nc_intrahem_lesioned.mat'),'-struct','S');
    case 'NC_getIntraNonlesionedHemPattern'   % get pattern of correlations for pairs within lesioned hemisphere
        % 0. Load data
        D       = load(fullfile(rootDir,'rs_preprocess.mat')); 
        D       = rs_imana_new('PP_removeInf',D);

        % 0. Obtain required data subset        
        %       - exclude unwanted ROIs
        inclROI = [1:5];        %no parietal or v1 rois         
        idx = ismember(D.regType(:,1),inclROI) & ismember(D.regType(:,2),inclROI);        
        D   = getrow(D,idx);        
        
        %       - get intrahemispheric fisherz correlations
        D   = getrow(D,D.hl~=0);
        Dc  = getrow(D,D.control==1);
        Dp  = getrow(D,D.control==0 & D.lesionSide==D.hl);
        D   = addstruct(Dp,Dc);
        disp(what);
        
        S   = [];
        % Loop over subj
        for sn=unique(D.subj_name)'
            Ds = getrow(D,strcmp(D.subj_name,sn));
            
            % Loop over week
            for w=unique(Ds.week)'
                Dw  = getrow(Ds,Ds.week==w);
                fprintf('%s - W%d\n',Dw.subj_name{1},Dw.week(1));
                
                M   = zeros(length(inclROI),length(inclROI));     % matrix of corr. for all regions
                for i=1:size(Dw.regType,1)
                    M(Dw.regType(i,1),Dw.regType(i,2)) = Dw.fz_r(i);
                end;
                M   = M+M';
                
                % add data to new structure
                Si      = getrow(Dw,1);
                Si      = rmfield(Si,{'rType','rpair','r'});
                Si.M    = squareform(M);
                S       = addstruct(S,Si);
            end;
        end;
        save(fullfile(rootDir,'nc_intrahem_nonlesioned.mat'),'-struct','S');         
        
    case 'NC_getInterHemPatternSH'      % get pattern of correlations for pairs between hemispheres (split half)
        % 0. Load data
        D       = load(fullfile(rootDir,'rs_preprocess_splithalf.mat')); 
        D       = rs_imana_new('PP_removeInf',D);

        % 0. Obtain required data subset        
        %       - exclude unwanted ROIs
        inclROI = [1:5];        %no parietal or v1 rois         
        idx = ismember(D.regType(:,1),inclROI) & ismember(D.regType(:,2),inclROI);        
        D   = getrow(D,idx);        
        
        %       - get interhemispheric fisherz correlations
        D   = getrow(D,D.hl==0);
        disp(what);
        
        S   = [];
        % Loop over subj
        for sn=unique(D.subj_name)'
            Ds = getrow(D,strcmp(D.subj_name,sn));
            
            % Loop over week
            for w=unique(Ds.week)'
                Dw  = getrow(Ds,Ds.week==w);
                fprintf('%s - W%d\n',Dw.subj_name{1},Dw.week(1));
                
                M = zeros(length(inclROI),length(inclROI),2);     % matrix of corr. for all regions
                
                for s=1:2
                    Dsp  = getrow(Dw,Dw.split==s);
                    
                    % get pattern for current split
                    for i=1:size(Dsp.regType,1)
                        M(Dsp.regType(i,1),Dsp.regType(i,2),s) = Dsp.fz_r(i);
                    end;
                end;
                    
                % get patterns for splits 1 and 2
                M1 = nonsym_squareform(M(:,:,1));
                M2 = nonsym_squareform(M(:,:,2));
                
                % add data to new structure
                Si          = getrow(Dw,1);
                Si          = rmfield(Si,{'rType','rpair','r'});
                Si.M1       = M1;
                Si.M2       = M2;
                Si.sh_r     = corr(M1',M2');
                Si.sh_fz_r  = fisherz(Si.sh_r);
                S           = addstruct(S,Si);
            end;
        end;
        save(fullfile(rootDir,'nc_interhem_splithalf.mat'),'-struct','S'); 
    case 'NC_getIntraHemLesionedSH'     % get pattern of correlations for pairs between hemispheres (split half)
        % 0. Load data
        D       = load(fullfile(rootDir,'rs_preprocess_splithalf.mat')); 
        D       = rs_imana_new('PP_removeInf',D);

        % 0. Obtain required data subset        
        %       - exclude unwanted ROIs
        inclROI = [1:5];        %no parietal or v1 rois         
        idx = ismember(D.regType(:,1),inclROI) & ismember(D.regType(:,2),inclROI);        
        D   = getrow(D,idx);        
        
        %       - get intrahemispheric fisherz correlations
        D   = getrow(D,D.hl~=0);
        Dc  = getrow(D,D.control==1);
        Dp  = getrow(D,D.control==0 & D.lesionSide==D.hl);
        D   = addstruct(Dp,Dc);
        disp(what);
        
        S   = [];
        % Loop over subj
        for sn=unique(D.subj_name)'
            Ds = getrow(D,strcmp(D.subj_name,sn));
            
            % Loop over week
            for w=unique(Ds.week)'
                Dw  = getrow(Ds,Ds.week==w);
                fprintf('%s - W%d\n',Dw.subj_name{1},Dw.week(1));
                
                M = zeros(length(inclROI),length(inclROI),2);     % matrix of corr. for all regions
                
                for s=1:2
                    Dsp  = getrow(Dw,Dw.split==s);
                    
                    % get pattern for current split
                    for i=1:size(Dsp.regType,1)
                        M(Dsp.regType(i,1),Dsp.regType(i,2),s) = Dsp.fz_r(i);
                    end;
                end;
                    
                % get patterns for splits 1 and 2
                M(:,:,1)    = M(:,:,1)+M(:,:,1)';
                M(:,:,2)    = M(:,:,2)+M(:,:,2)';
                M1          = squareform(M(:,:,1));
                M2          = squareform(M(:,:,2));
                
                % add data to new structure
                Si          = getrow(Dw,1);
                Si          = rmfield(Si,{'rType','rpair','r'});
                Si.M1       = M1;
                Si.M2       = M2;
                Si.sh_r     = corr(M1',M2');
                Si.sh_fz_r  = fisherz(Si.sh_r);
                S           = addstruct(S,Si);
            end;
        end;
        save(fullfile(rootDir,'nc_intrahem_lesioned_splithalf.mat'),'-struct','S');   
    case 'NC_getIntraHemNonlesionedSH'  % get pattern of correlations for pairs between hemispheres (split half)
        % 0. Load data
        D       = load(fullfile(rootDir,'rs_preprocess_splithalf.mat')); 
        D       = rs_imana_new('PP_removeInf',D);

        % 0. Obtain required data subset        
        %       - exclude unwanted ROIs
        inclROI = [1:5];        %no parietal or v1 rois         
        idx = ismember(D.regType(:,1),inclROI) & ismember(D.regType(:,2),inclROI);        
        D   = getrow(D,idx);        
        
        %       - get intrahemispheric fisherz correlations
        D   = getrow(D,D.hl~=0);
        Dc  = getrow(D,D.control==1);
        Dp  = getrow(D,D.control==0 & D.lesionSide~=D.hl);
        D   = addstruct(Dp,Dc);
        disp(what);
        
        S   = [];
        % Loop over subj
        for sn=unique(D.subj_name)'
            Ds = getrow(D,strcmp(D.subj_name,sn));
            
            % Loop over week
            for w=unique(Ds.week)'
                Dw  = getrow(Ds,Ds.week==w);
                fprintf('%s - W%d\n',Dw.subj_name{1},Dw.week(1));
                
                M = zeros(length(inclROI),length(inclROI),2);     % matrix of corr. for all regions
                
                for s=1:2
                    Dsp  = getrow(Dw,Dw.split==s);
                    
                    % get pattern for current split
                    for i=1:size(Dsp.regType,1)
                        M(Dsp.regType(i,1),Dsp.regType(i,2),s) = Dsp.fz_r(i);
                    end;
                end;
                    
                % get patterns for splits 1 and 2
                M(:,:,1)    = M(:,:,1)+M(:,:,1)';
                M(:,:,2)    = M(:,:,2)+M(:,:,2)';
                M1          = squareform(M(:,:,1));
                M2          = squareform(M(:,:,2));
                
                % add data to new structure
                Si          = getrow(Dw,1);
                Si          = rmfield(Si,{'rType','rpair','r'});
                Si.M1       = M1;
                Si.M2       = M2;
                Si.sh_r     = corr(M1',M2');
                Si.sh_fz_r  = fisherz(Si.sh_r);
                S           = addstruct(S,Si);
            end;
        end;
        save(fullfile(rootDir,'nc_intrahem_nonlesioned_splithalf.mat'),'-struct','S');           
        
    case 'REL_withinSubj'               % within subject correlations for patterns over weeks
        % 0. Get data
        %   - data could be intra, inter lesionsed, inter nonlesioned
        D       = varargin{1};
        
        %   - mark weeks 1-4 as early, 24-52 as late
        week    = {[1,4],[24,52]};
        D.stage = zeros(length(D.subj_name),1);
        for i=1:2
            D.stage(ismember(D.week,week{i})) = i;
        end;
        D = getrow(D,D.stage~=0);
        
        %   - get only those subjs that have measurements for early and
        %   late
        [x,sn]      = pivottable(D.subj_name,D.stage,D.stage,'length(unique(x))');
        x(isnan(x)) = 0;
        x           = all(x,2);
        D           = getrow(D,ismember(D.subj_name,sn(x)));
        D = tapply(D,{'subj_name','stage','control'},{'M','nanmean(x,1)'});
        
        %   - calculate correlations between early and late stages
        S = [];
        for sn=unique(D.subj_name)'
            Ds1 = getrow(D,ismember(D.subj_name,sn));
            Ds2 = getrow(D,~ismember(D.subj_name,sn) & D.control==mean(Ds1.control));
            
            M1 = Ds1.M(Ds1.stage==1,:);   % early
            M2 = Ds1.M(Ds1.stage==2,:);   % late
            
            rEarly  = mean(fisherz(corr(M1',Ds2.M(Ds2.stage==1,:)')));
            rLate   = mean(fisherz(corr(M2',Ds2.M(Ds2.stage==2,:)')));

            Si          = getrow(Ds1,1);
            Si          = rmfield(Si,{'stage','M'});
            Si.r        = corr(M1',M2');
            Si.fz_r     = fisherz(Si.r);
            Si.fz_early = rEarly;
            Si.fz_late  = rLate;
            S           = addstruct(S,Si);
        end;
        varargout = {S};
    case 'REL_withinAll'                % calculate within subject correlations (early/late after stroke)
        fName = {'nc_interhem','nc_intrahem_lesioned','nc_intrahem_nonlesioned'};
        
        S = [];
        for i=1:length(fName)
            D = load(fullfile(rootDir,sprintf('%s.mat',fName{i})));
            D = rs_imana_new('REL_withinSubj',D);
            
            Si      = D;
            Si.type = zeros(length(D.subj_name),1)+i;
            S       = addstruct(S,Si);
            
            disp(fName{i});
            x = pivottable(D.subj_name,D.control,D.fz_r,'nanmean');
            ttest(x(:,1),x(:,2),2,'independent');
        end;
        barplot(S.type,S.fz_r,'split',S.control);
        
        [~,~,sn] = unique(S.subj_name);
        anovaMixed(S.fz_r,sn,'between',S.control,{'grp'},'within',S.type,{'type'});
    case 'REL_fingerprint'              % check to see if a finger print exists for the participants
        fName = {'nc_interhem','nc_intrahem_lesioned','nc_intrahem_nonlesioned'};
        
        S = [];
        for i=1:length(fName)
            D = load(fullfile(rootDir,sprintf('%s.mat',fName{i})));
            D = rs_imana_new('REL_withinSubj',D);
            
            Si      = D;
            Si.type = zeros(length(D.subj_name),1)+i;
            S       = addstruct(S,Si);
        end;
        
        T = [];
        l = {'fz_r','fz_early','fz_late'};
        for i=1:length(l)
            Ti      = S;
            Ti.fp   = zeros(length(S.subj_name),1)+i;
            Ti.fz   = Ti.(l{i});
            T       = addstruct(T,Ti);
        end;
        
        subplot(121);
        barplot(T.type,T.fz,'split',T.fp,'subset',T.control==1);
        title('control');
        subplot(122);
        barplot(T.type,T.fz,'split',T.fp,'subset',T.control==0);
        title('patient');
        
        [~,~,T.sn] = unique(T.subj_name);
        T1 = getrow(T,T.control==1);
        disp('controls');
        anovaMixed(T1.fz,T1.sn,'within',[T1.type T1.fp],{'fingerprint','type'});
        
        T2 = getrow(T,T.control==0);
        disp('patients');
        anovaMixed(T2.fz,T2.sn,'within',[T2.type T2.fp],{'fingerprint','type'});
    case 'PP_getSubset'                 % get subset of data from the main data structure
        exclROI     = [6:8];                % regions to exclude (no parietal or V1)
        exclSubj    = {''};                 % remove these participants
        type        = 'lesioned';           % lesioned/nonlesioned/inter/all
        fname       = 'rs_preprocess.mat';  % over-ride file name to open a different structure
        lesionType  = [];                   % all lesion types
        vararginoptions(varargin,{'exclROI','exclSubj','type','fname','lesionType'});
        
        % 1. Load data for split-half reliabilities
        %       - only include required ROIs
        D   = load(fullfile(ppDir,fname));
        idx = ~(ismember(D.regType(:,1),exclROI) | ismember(D.regType(:,2),exclROI));
        D   = getrow(D,idx);
        
        %       - remove these subjects
        idx = ~ismember(D.subj_name,exclSubj);
        D   = getrow(D,idx);

        %       - get rows corresponding only to the type
        switch(type)
            case 'lesioned'
                Dc  = getrow(D,D.control==1 & D.hl~=0);                 % controls, intra-cortical only
                Dp  = getrow(D,D.control==0 & D.lesionSide==D.hl);      % patients, lesioned hem, intra-cortical only.         
            case 'nonlesioned'
                Dc  = getrow(D,D.control==1 & D.hl~=0);                 % controls, intra-cortical only
                Dp  = getrow(D,D.control==0 & D.lesionSide~=D.hl);      % patients, nonlesioned hem, intra-cortical only.         
            case 'inter'
                Dc  = getrow(D,D.control==1 & D.hl==0);                 % controls, inter-cortical only
                Dp  = getrow(D,D.control==0 & D.hl==0);                 % patients, inter-cortical only                 
            case 'all'
                Dc  = getrow(D,D.control==1);                           % all controls
                Dp  = getrow(D,D.control==0);                           % all patients
        end;
        
        %       - only get specified lesion types
        if ~isempty(lesionType)
            Dp  = getrow(Dp,ismember(Dp.lesionType,lesionType));
        end;
        
        % 2. return structure
        T   = addstruct(Dc,Dp);        
        varargout = {T};
    case 'PP_patientSplit'              % get subset split of patients
%         Dp = varargin{1};
%         vararginoptions({varargin{2:end}},{'lesionType'});
%         
%         % 0. Avg controls (all) and patients (per week and according to severity)
%         C           = tapply(D,{'subj_name'},{'mat','mean(x,1)','subset',D.control==1});
% 
%         D.mvcNorm   = -D.mvcNorm;                       % for some reason mvcNorm is negative (to be fixed)
%         [x,sn]      = pivottable(D.subj_name,D.week,D.mvcNorm,'mean(x)','subset',D.control==0);
%         x           = nanmean(x(:,1:2),2);

    case 'FIG_interhemSomat'            % somatotopic ordering of correlations between hemispheres    
        S = load(fullfile(rootDir,'nc_interhem.mat'));
        
        hom = 1:6:25;                       % homologous interhemispheric areas
        het = find(~ismember(1:25,hom));    % heterologous interhemispheric areas
        
        % get required roi pairs
        S.hom = nanmean(S.M(:,hom),2);
        S.het = nanmean(S.M(:,het),2);
        
        % make plot for hom/het pairs
        figure;
        subplot(121);
        lineplot(S.week,[S.hom S.het],'subset',S.control==1,'style_thickline','leg',{'hom','het'});
        title('controls');
        ylim([0.2 1.1]);
        subplot(122);
        lineplot(S.week,[S.hom S.het],'subset',S.control==0,'style_thickline','leg',{'hom','het'});
        title('patients');
        ylim([0.2 1.1]);
        
        % do ttest for controls (hom vs het)
        x=pivottablerow(S.subj_name,[S.hom S.het],'nanmean(x,1)','subset',S.control==1);
        ttest(x(:,1),x(:,2),2,'paired');        
        
        % save for statistics
        type    = zeros(length(S.control),1);
        D1      = S;
        D1.type = type + 1;
        D1.fz   = D1.hom;
        D2      = S;
        D2.type = type + 2;
        D2.fz   = D2.het;
        D       = addstruct(D1,D2);
        D = tapply(D,{'subj_name','week','control','type'},{'fz','nanmean'});
        dsave(fullfile(rDir,[what,'.dat']),D);
    case 'FIG_inter_splithalf_reliability'    % split half reliability of patterns within each individual/session
        S = load(fullfile(rootDir,'nc_interhem_splithalf.mat'));        
        lineplot(S.week,S.sh_fz_r,'split',S.control,'plotfcn','fisherinv(nanmean(x))','style_thickline','leg','auto');
        
        % save for statistics
        D = tapply(S,{'subj_name','week','control'},{'sh_fz_r','nanmean'});
        dsave(fullfile(rDir,[what,'.dat']),D);
    case 'FIG_intraLesioned_splithalf_reliability'    % split half reliability of patterns within each individual/session
        S = load(fullfile(rootDir,'nc_intraHem_lesioned_splithalf.mat'));        
        lineplot(S.week,S.sh_fz_r,'split',S.control,'plotfcn','fisherinv(nanmean(x))','style_thickline','leg','auto');
        
        % save for statistics
        D = tapply(S,{'subj_name','week','control'},{'sh_fz_r','nanmean'});
        dsave(fullfile(rDir,[what,'.dat']),D);        
    case 'FIG_intraNonlesioned_splithalf_reliability'    % split half reliability of patterns within each individual/session
        S = load(fullfile(rootDir,'nc_intraHem_nonlesioned_splithalf.mat'));        
        lineplot(S.week,S.sh_fz_r,'split',S.control,'plotfcn','fisherinv(nanmean(x))','style_thickline','leg','auto');
        
        % save for statistics
        D = tapply(S,{'subj_name','week','control'},{'sh_fz_r','nanmean'});
        dsave(fullfile(rDir,[what,'.dat']),D);                
        
    case 'NC_getIntra'          % get intra hemispheric correlations for controls and patients
        type = 'ens';
        vararginoptions(varargin,{'D','lesionSide',});
        
        D       = load(fullfile(rootDir,'rs_preprocess.mat')); 
        exclROI = [6 7 8];        %no parietal or v1 rois 
        
        % exclude ROIs
        idx = ~(ismember(D.regType(:,1),exclROI) | ismember(D.regType(:,2),exclROI));
        D   = getrow(D,idx);        
        Dp  = getrow(D,D.control==0 & D.lesionSide==D.hl);      % get patients, lesioned hem, intra-cortical only. 
        Dc  = getrow(D,D.control==1 & D.hl~=0);                 % get controls, intra-cortical only
        D   = addstruct(Dc,Dp);
        
        S   = [];
        for sn=unique(D.subj_name)'
            Ds = getrow(D,strcmp(D.subj_name,sn));
            for w=unique(Ds.week)'
                Dw  = getrow(Ds,Ds.week==w);
                fprintf('%s - W%d\n',Dw.subj_name{1},Dw.week(1));
                
                M   = zeros(8,8);     % matrix of corr. for all regions
                for i=1:size(Dw.regType,1)
                    M(Dw.regType(i,1),Dw.regType(i,2)) = Dw.fz_r(i);
                end;
                M = M+M';       % symmetricizing matrix
                idx = find(~ismember(1:8,exclROI));
                M = M(idx,idx);
                
                % add data to new structure
                Si      = getrow(Dw,1);
                Si      = rmfield(Si,{'rType','rpair','r'});
                Si.M    = nonsym_squareform(M);
                S       = addstruct(S,Si);
            end;
        end;
        save(fullfile(rootDir,'rs_networkCorr_lesionedIntra.mat'),'-struct','S');         
        
    case 'cf_connectivity'
    D   = load(fullfile(rootDir,'rs_preprocess.mat'));  
    
    % only take subcorticals for now
      D   = getrow(D,(D.lesionType~=2)); %no cortical
      D   = getrow(D,(D.lesionType~=4)); %no mixed
    
%     or define threshold
%     D   = getrow(D,D.lesion_PrCG<=0.05);
    
%     only take chosen ones
     %subj = {'CU_2925','JHU_2531','UZ_2654','UZ_3239','UZ_3240','UZ_3241',...
     %   'UZ_3243','UZ_3246','UZ_3247','UZ_3248','CU_1002','JHP_1001',...
     %   'JHP_1002','JHP_1004','UZP_1001','UZP_1002','UZP_1004','UZP_1005',...
     %   'UZP_1006','UZP_1007','UZP_1008'};
     %D = getrow(D,ismember(D.subj_name,subj));
    
%     only exclude the bad ones 
%    subj = {'CU_2697','JHU_2374','UZ_3228'};
%       D = getrow(D,~ismember(D.subj_name,subj));

    
%  Roi pairs
     M1_M1    = getrow(D,D.rType==23);        %[2 10] 
     Pmd_Pmd  = getrow(D,D.rType==37);        %[3 11]
     Pmv_Pmv  = getrow(D,D.rType==50);        %[4 12]
     SMA_SMA  = getrow(D,D.rType==62);        %[5 13]
     V1_V1    = getrow(D,D.rType==73);        %[6 14]

     
    s   = style_sheet(sty_2grp_ac,'leg',{'patients','control'},'leglocation','northoutside');
     
    %interhemispheric connectivity 
    h = make_figure;
    subplot(3,2,1)
    title('M1-M1 subcortical')
    %lineplot(M1_M1.week,M1_M1.fz_r,'plotfcn','nanmean','errorfcn','stderr','split',M1_M1.control,'CAT',s(1).CAT,s(1).PLOT{:});
    lineplot(M1_M1.week,M1_M1.r,'plotfcn','nanmean','errorfcn','stderr','split',M1_M1.control,'CAT',s(1).CAT,s(1).PLOT{:});
    ylabel('correlation index')
    xlabel('weeks')
    subplot(3,2,2)
    title('Pmd-Pmd subcortical')
    %lineplot(Pmd_Pmd.week,Pmd_Pmd.fz_r,'plotfcn','nanmean','errorfcn','stderr','split',Pmd_Pmd.control,'CAT',s(1).CAT,s(1).PLOT{:});
    lineplot(Pmd_Pmd.week,Pmd_Pmd.r,'plotfcn','nanmean','errorfcn','stderr','split',Pmd_Pmd.control,'CAT',s(1).CAT,s(1).PLOT{:});
    ylabel('correlation index')
    xlabel('weeks')
    subplot(3,2,3)
    title('Pmv-Pmv subcortical')
    %lineplot(Pmv_Pmv.week,Pmv_Pmv.fz_r,'plotfcn','nanmean','errorfcn','stderr','split',Pmv_Pmv.control,'CAT',s(1).CAT,s(1).PLOT{:});
    lineplot(Pmv_Pmv.week,Pmv_Pmv.r,'plotfcn','nanmean','errorfcn','stderr','split',Pmv_Pmv.control,'CAT',s(1).CAT,s(1).PLOT{:});
    ylabel('correlation index')
    xlabel('weeks')
    subplot(3,2,4)
    title('SMA-SMA subcortical')
    %lineplot(SMA_SMA.week,SMA_SMA.fz_r,'plotfcn','nanmean','errorfcn','stderr','split',SMA_SMA.control,'CAT',s(1).CAT,s(1).PLOT{:});
    lineplot(SMA_SMA.week,SMA_SMA.r,'plotfcn','nanmean','errorfcn','stderr','split',SMA_SMA.control,'CAT',s(1).CAT,s(1).PLOT{:});
    ylabel('correlation index')
    xlabel('weeks') 
    subplot(3,2,5)
    title('V1-V1 subcortical')
    %lineplot(V1_V1.week,V1_V1.fz_r,'plotfcn','nanmean','errorfcn','stderr','split',V1_V1.control,'CAT',s(1).CAT,s(1).PLOT{:});
    lineplot(V1_V1.week,V1_V1.r,'plotfcn','nanmean','errorfcn','stderr','split',V1_V1.control,'CAT',s(1).CAT,s(1).PLOT{:});
    ylabel('correlation index')
    xlabel('weeks') 
    
    %% interhemispheric individual patients
    p   = style_sheet(sty_multiple,'leg','auto','leglocation','northoutside');
    
    h = make_figure;
    subplot(3,2,1)
    title('M1-M1 subcortical')
    %lineplot(M1_M1.week,M1_M1.fz_r,'plotfcn','nanmean','errorfcn','stderr','split',M1_M1.subj_name,'subset',M1_M1.control==0);
    lineplot(M1_M1.week,M1_M1.r,'plotfcn','nanmean','errorfcn','stderr','split',M1_M1.subj_name,'subset',M1_M1.control==0);
    ylabel('correlation index')
    xlabel('weeks')
    subplot(3,2,2)
    title('Pmd-Pmd subcortical')
    lineplot(Pmd_Pmd.week,Pmd_Pmd.fz_r,'plotfcn','nanmean','errorfcn','stderr','split',Pmd_Pmd.subj_name,'subset',Pmd_Pmd.control==0);
    %lineplot(Pmd_Pmd.week,Pmd_Pmd.r,'plotfcn','nanmean','errorfcn','stderr','split',Pmd_Pmd.subj_name,'subset',Pmd_Pmd.control==0);
    ylabel('correlation index')
    xlabel('weeks')
    subplot(3,2,3)
    title('Pmv-Pmv subcortical')
    %lineplot(Pmv_Pmv.week,Pmv_Pmv.fz_r,'plotfcn','nanmean','errorfcn','stderr','split',Pmv_Pmv.subj_name,'subset',Pmv_Pmv.control==0);
    lineplot(Pmv_Pmv.week,Pmv_Pmv.r,'plotfcn','nanmean','errorfcn','stderr','split',Pmv_Pmv.subj_name,'subset',Pmv_Pmv.control==0);
    ylabel('correlation index')
    xlabel('weeks')
    subplot(3,2,4)
    title('SMA-SMA subcortical')
    %lineplot(SMA_SMA.week,SMA_SMA.fz_r,'plotfcn','nanmean','errorfcn','stderr','split',SMA_SMA.subj_name,'subset',SMA_SMA.control==0);
    lineplot(SMA_SMA.week,SMA_SMA.r,'plotfcn','nanmean','errorfcn','stderr','split',SMA_SMA.subj_name,'subset',SMA_SMA.control==0);
    ylabel('correlation index')
    xlabel('weeks')
    subplot(3,2,5)
    title('V1-V1 subcortical')
    %lineplot(V1_V1.week,V1_V1.fz_r,'plotfcn','nanmean','errorfcn','stderr','split',V1_V1.subj_name,'subset',V1_V1.control==0);
    lineplot(V1_V1.week,V1_V1.r,'plotfcn','nanmean','errorfcn','stderr','split',V1_V1.subj_name,'subset',V1_V1.control==0);
    ylabel('correlation index')
    xlabel('weeks')
    
%% intracortical connectivity/ lesioned/dominant
    
     h = make_figure;
     subplot(1,4,1)
     title('M1-SMA subcortical')
     M1_SMA = getrow(T, T.regType(:,1)==2 & T.regType(:,2)==5);
     lineplot(M1_SMA.week,M1_SMA.r,'plotfcn','nanmean','errorfcn','stderr','split',M1_SMA.control,'CAT',s(1).CAT,s(1).PLOT{:}); 
     %individuals
     %lineplot(M1_SMA.week,M1_SMA.r,'plotfcn','nanmean','errorfcn','stderr','split',M1_SMA.subj_name,'subset',M1_SMA.control==0); 
     ylabel('correlation index')
     xlabel('weeks')
     subplot(1,4,2)
     title('M1-Pmd subcortical')
     M1_Pmd = getrow(T, T.regType(:,1)==2 & T.regType(:,2)==3);
     lineplot(M1_Pmd.week,M1_Pmd.r,'plotfcn','nanmean','errorfcn','stderr','split',M1_Pmd.control,'CAT',s(1).CAT,s(1).PLOT{:});
     %individuals
     %lineplot(M1_Pmd.week,M1_Pmd.r,'plotfcn','nanmean','errorfcn','stderr','split',M1_Pmd.subj_name,'subset',M1_Pmd.control==0);
     ylabel('correlation index')
     xlabel('weeks')
     subplot(1,4,3)
     title('M1-Pmv subcortical')
     M1_Pmv = getrow(T, T.regType(:,1)==2 & T.regType(:,2)==4);
     lineplot(M1_Pmv.week,M1_Pmv.r,'plotfcn','nanmean','errorfcn','stderr','split',M1_Pmv.control,'CAT',s(1).CAT,s(1).PLOT{:});
     %individuals
     %lineplot(M1_Pmv.week,M1_Pmv.r,'plotfcn','nanmean','errorfcn','stderr','split',M1_Pmv.subj_name,'subset',M1_Pmd.control==0); 
     ylabel('correlation index')
     xlabel('weeks')
     subplot(1,4,4)
     title('M1-S1 subcortical')
     M1_S1 = getrow(T, T.regType(:,1)==1 & T.regType(:,2)==2);
     lineplot(M1_S1.week,M1_S1.r,'plotfcn','nanmean','errorfcn','stderr','split',M1_S1.control,'CAT',s(1).CAT,s(1).PLOT{:});
     %individuals
     %lineplot(M1_S1.week,M1_S1.r,'plotfcn','nanmean','errorfcn','stderr','split',M1_S1.subj_name,'subset',M1_S1.control==0);
     ylabel('correlation index')
     xlabel('weeks')
     
%% intracortical connectivity/ NON lesioned/dominant
     T = getrow(D, D.lesionSide==D.hl);     
     h = make_figure;
     subplot(1,4,1)
     title('M1-SMA subcortical')
     M1_SMA = getrow(T, T.regType(:,1)==2 & T.regType(:,2)==5);
     lineplot(M1_SMA.week,M1_SMA.r,'plotfcn','nanmean','errorfcn','stderr','split',M1_SMA.control); ylabel('correlation index')
     xlabel('weeks')
     subplot(1,4,2)
     title('M1-Pmd subcortical')
     M1_Pmd = getrow(T, T.regType(:,1)==2 & T.regType(:,2)==3);
     lineplot(M1_Pmd.week,M1_Pmd.r,'plotfcn','nanmean','errorfcn','stderr','split',M1_Pmd.control,'CAT',s(1).CAT,s(1).PLOT{:});
     ylabel('correlation index')
     xlabel('weeks')
     subplot(1,4,3)
     title('M1-Pmv subcortical')
     M1_Pmv = getrow(T, T.regType(:,1)==2 & T.regType(:,2)==4);
     lineplot(Pmv.week,Pmv.r,'plotfcn','nanmean','errorfcn','stderr','split',M1.control,'CAT',s(1).CAT,s(1).PLOT{:});
     ylabel('correlation index')
     xlabel('weeks')
     subplot(1,4,4)
     title('M1-S1 subcortical')
     M1_S1 = getrow(T, T.regType(:,1)==2 & T.regType(:,2)==1);
     lineplot(M1_S1.week,M1_S1.r,'plotfcn','nanmean','errorfcn','stderr','split',M1_S1.control,'CAT',s(1).CAT,s(1).PLOT{:});
     ylabel('correlation index')
     xlabel('weeks')

keyboard;       
    case 'cf_conn_behav_inter'
    D   = load(fullfile(rootDir,'rs_preprocess.mat'));  
    
    % only take subcorticals for now
    D   = getrow(D,(D.lesionType~=2)); %no cortical
    D   = getrow(D,(D.lesionType~=4)); %no mixed
    
    M1_M1    = getrow(D,D.rType==23);        %[2 10] 
    
    s   = style_sheet(sty_2grp_ac,'leg',{'patients','control'},'leglocation','northoutside');
    
    h = make_figure;
    for w = 1:length(unique(M1_M1.week))
    we = M1_M1.week(w);
    subplot(1,5,w)
    scatterplot(M1_M1.fz_r,M1_M1.mvcNorm,'regression','linear','printcorr','split',M1_M1.control,'subset',M1_M1.week==we,'CAT',s(1).CAT,s(1).PLOT{:});
    xlabel('M1 M1 correlation')
    ylabel('mvcNorm')
    title(we)
    end   
    
    h = make_figure;
    for w = 1:length(unique(M1_M1.week))
    we = M1_M1.week(w);
    subplot(1,5,w)
    scatterplot(M1_M1.fz_r,M1_M1.ensOverall,'regression','linear','printcorr','split',M1_M1.control,'subset',M1_M1.week==we,'CAT',s(1).CAT,s(1).PLOT{:});
    xlabel('M1_M1 correlation')
    ylabel('enslaving')
    title(we)
    end
    
    h = make_figure;
    for w = 1:length(unique(M1_M1.week))
    we = M1_M1.week(w);
    subplot(1,5,w)
    scatterplot(M1_M1.fz_r,M1_M1.mmOverall,'regression','linear','printcorr','split',M1_M1.control,'subset',M1_M1.week==we,'CAT',s(1).CAT,s(1).PLOT{:});
    xlabel('M1_M1 correlation')
    ylabel('mirror movements')
    title(we)
    end
    
    h = make_figure;
    for w = 1:length(unique(M1_M1.week))
    we = M1_M1.week(w);
    subplot(1,5,w)
    scatterplot(M1_M1.fz_r,M1_M1.FM,'regression','linear','printcorr','subset',M1_M1.week==we & M1_M1.control==0,'CAT',s(1).CAT,s(1).PLOT{:});
    xlabel('M1_M1 correlation')
    ylabel('FM')
    title(we)
    end
    
    h = make_figure;
    corr_w1 = getrow(M1_M1,M1_M1.week==1); 
    behav_52 = getrow(M1_M1,M1_M1.week==52);
    subplot(1,4,1)
    scatterplot(corr_w1.fz_r,behav_52.mvcNorm,'regression','linear','printcorr','subset',corr_w1.control==0,'CAT',s(1).CAT,s(1).PLOT{:});
    xlabel('M1_M1 acute')
    ylabel('mvc')
    subplot(1,4,2)
    scatterplot(corr_w1.fz_r,behav_52.ensOverall,'regression','linear','printcorr','subset',corr_w1.control==0,'CAT',s(1).CAT,s(1).PLOT{:});
    xlabel('M1_M1 acute')
    ylabel('enslaving')
    subplot(1,4,3)
    scatterplot(corr_w1.fz_r,behav_52.mmOverall,'regression','linear','printcorr','subset',corr_w1.control==0,'CAT',s(1).CAT,s(1).PLOT{:});
    xlabel('M1_M1 acute')
    ylabel('mirror movements')
    subplot(1,4,4)
    scatterplot(corr_w1.fz_r,behav_52.FM,'regression','linear','printcorr','subset',corr_w1.control==0,'CAT',s(1).CAT,s(1).PLOT{:});
    xlabel('M1_M1 acute')
    ylabel('FM chronic')
    
    
    keyboard;
    case 'cf_conn_M1_CST'   
    D   = load(fullfile(rootDir,'rs_preprocess.mat'));  
    
    % only take subcorticals for now
    D   = getrow(D,(D.lesionType~=2)); %no cortical
    D   = getrow(D,(D.lesionType~=4)); %no mixed
    
    M1_M1    = getrow(D,D.rType==23);        %[2 10] 
     
    h = make_figure;
    M1_M1   = getrow(M1_M1,M1_M1.control==0);
    for w = 1:length(unique(M1_M1.week))
    we = M1_M1.week(w);
    subplot(1,5,w)
    scatterplot(M1_M1.lesion_CST_total,M1_M1.r,'regression','linear','printcorr','subset',M1_M1.week==we,'markercolor',sty_2grp_ac,'markerfill',sty_2grp_ac);
    end
    
    h = make_figure;
    M1_M1   = getrow(M1_M1,M1_M1.control==0);
    M1_M1   = getrow(M1_M1,M1_M1.lesion_CST_total>0);
    for w = 1:length(unique(M1_M1.week))
    we = M1_M1.week(w);
    subplot(1,5,w)
    scatterplot(M1_M1.lesion_CST_total,M1_M1.r,'regression','linear','printcorr','subset',M1_M1.week==we,'markercolor',sty_2grp_ac,'markerfill',sty_2grp_ac);
    end    
    case 'make_networkCorr'         % correlate motor networks for each subj/week
        D       = load(fullfile(rootDir,'rs_preprocess.mat')); 
        exclROI = [6 7 8];        %no parietal or v1 rois 
        
        % exclude ROIs
        idx = ~(ismember(D.regType(:,1),exclROI) | ismember(D.regType(:,2),exclROI));
        D   = getrow(D,idx);        
        D  = getrow(D,D.hl==0);      % inter-cortical only. 
        
        S   = [];
        for sn=unique(D.subj_name)'
            Ds = getrow(D,strcmp(D.subj_name,sn));
            for w=unique(Ds.week)'
                Dw  = getrow(Ds,Ds.week==w);
                fprintf('%s - W%d\n',Dw.subj_name{1},Dw.week(1));
                
                M   = zeros(10,10);     % matrix of corr. for all regions
                for i=1:size(Dw.rpair,1)
                    M(Dw.rpair(i,1),Dw.rpair(i,2)-3) = Dw.fz_r(i);
                end;
                M = M+M';       % symmetricizing matrix
                
                % add data to new structure
                Si      = getrow(Dw,1);
                Si      = rmfield(Si,{'rType','rpair','r'});
                Si.M    = nonsym_squareform(M);
                Si.mat  = nonsym_squareform(M(1:5,6:10));
                S       = addstruct(S,Si);
            end;
        end;
        save(fullfile(rootDir,'rs_networkCorr.mat'),'-struct','S');   
        
    case 'make_networkCorr_lesionedIntra'         % correlate motor networks, only lesioned hem in patients
        D       = load(fullfile(rootDir,'rs_preprocess.mat')); 
        exclROI = [6 7 8];        %no parietal or v1 rois 
        
        % exclude ROIs
        idx = ~(ismember(D.regType(:,1),exclROI) | ismember(D.regType(:,2),exclROI));
        D   = getrow(D,idx);        
        Dp  = getrow(D,D.control==0 & D.lesionSide==D.hl);      % get patients, lesioned hem, intra-cortical only. 
        Dc  = getrow(D,D.control==1 & D.hl~=0);                 % get controls, intra-cortical only
        D   = addstruct(Dc,Dp);
        
        S   = [];
        for sn=unique(D.subj_name)'
            Ds = getrow(D,strcmp(D.subj_name,sn));
            for w=unique(Ds.week)'
                Dw  = getrow(Ds,Ds.week==w);
                fprintf('%s - W%d\n',Dw.subj_name{1},Dw.week(1));
                
                M   = zeros(8,8);     % matrix of corr. for all regions
                for i=1:size(Dw.regType,1)
                    M(Dw.regType(i,1),Dw.regType(i,2)) = Dw.fz_r(i);
                end;
                M = M+M';       % symmetricizing matrix
                idx = find(~ismember(1:8,exclROI));
                M = M(idx,idx);
                
                % add data to new structure
                Si      = getrow(Dw,1);
                Si      = rmfield(Si,{'rType','rpair','r'});
                Si.M    = nonsym_squareform(M);
                S       = addstruct(S,Si);
            end;
        end;
        save(fullfile(rootDir,'rs_networkCorr_lesionedIntra.mat'),'-struct','S');                        
    case 'make_networkCorr_nonlesionedIntra'         % correlate motor networks, only lesioned hem in patients
        D       = load(fullfile(rootDir,'rs_preprocess.mat')); 
        exclROI = [6 7 8];        %no parietal or v1 rois 
        
        % exclude ROIs
        idx = ~(ismember(D.regType(:,1),exclROI) | ismember(D.regType(:,2),exclROI));
        D   = getrow(D,idx);        
        Dp  = getrow(D,D.control==0 & D.hl~=0 & D.lesionSide~=D.hl);      % get patients, non-lesioned hem, intra-cortical only. 
        Dc  = getrow(D,D.control==1 & D.hl~=0);                 % get controls, intra-cortical only
        D   = addstruct(Dc,Dp);
        
        S   = [];
        for sn=unique(D.subj_name)'
            Ds = getrow(D,strcmp(D.subj_name,sn));
            for w=unique(Ds.week)'
                Dw  = getrow(Ds,Ds.week==w);
                fprintf('%s - W%d\n',Dw.subj_name{1},Dw.week(1));
                
                M   = zeros(8,8);     % matrix of corr. for all regions
                for i=1:size(Dw.regType,1)
                    M(Dw.regType(i,1),Dw.regType(i,2)) = Dw.fz_r(i);
                end;
                M = M+M';       % symmetricizing matrix
                idx = find(~ismember(1:8,exclROI));
                M = M(idx,idx);
                
                % add data to new structure
                Si      = getrow(Dw,1);
                Si      = rmfield(Si,{'rType','rpair','r'});
                Si.M    = nonsym_squareform(M);
                S       = addstruct(S,Si);
            end;
        end;
        save(fullfile(rootDir,'rs_networkCorr_nonlesionedIntra.mat'),'-struct','S');                                
        
    case 'plot_networkCorr_all'     % plot correlation across controls and patients
        D   = load(fullfile(rootDir,'rs_networkCorr.mat'));                
        %take only subgroup
         D   = getrow(D,(D.lesionType~=2)); %no cortical
         D   = getrow(D,(D.lesionType~=4)); %no mixed
         %D   = getrow(D,D.lesion_PrCG<=0.05);
        %sn  = {'UZ_3030','UZ_3226','UZ_3227'};
        %D   = getrow(D,~ismember(D.subj_name,sn));
        %D   = getrow(D,~(ismember(D.subj_name,'JHU_2713') & D.week==12));

        
        Dc  = getrow(D,D.control==1);
        Dp  = getrow(D,D.control~=1);
        
        % get avg matrices for each group
        Mc  = nanmean(Dc.M,1);
        Mc  = nonsym_squareform(Mc);
        
        
        DpW1 = getrow (Dp,Dp.week==1);
        Mp1 = nanmean(DpW1.M,1);
        Mp1 = nonsym_squareform(Mp1);
        DpW4 = getrow (Dp,Dp.week==4);
        Mp4 = nanmean(DpW4.M,1);
        Mp4 = nonsym_squareform(Mp4);
        DpW12 = getrow (Dp,Dp.week==12);
        Mp12 = nanmean(DpW12.M,1);
        Mp12 = nonsym_squareform(Mp12);
        DpW24 = getrow (Dp,Dp.week==24);
        Mp24 = nanmean(DpW24.M,1);
        Mp24 = nonsym_squareform(Mp24);
        DpW52 = getrow (Dp,Dp.week==52);
        Mp52 = nanmean(DpW52.M,1);
        Mp52 = nonsym_squareform(Mp52);
        
%         Mp = [];
%         for w=unique(Dp.week)'
%             Dw  = getrow(Dp,Dp.week==w);
%             Mpw  = nanmean(Dw.M,1);
%             Mpw  = nonsym_squareform(Mpw);
%             Mp = addstruct(Mp,Mpw);
%         end
        scale = [0 0.8];
        h=make_figure;
        subplot(231);
        colormap(hot); 
        imagesc(Mc,scale); colorbar;
        title('Controls');
        set_graphics(h,'xtick',1:16,'ytick',1:16, 'ax','square')
        subplot(232);
        colormap(hot); 
        imagesc(Mp1,scale); colorbar;
        title('Patients Week1');
        set_graphics(h,'xtick',1:16,'ytick',1:16, 'ax','square')
        subplot(233);
        colormap(hot); 
        imagesc(Mp4,scale); colorbar;
        title('Patients Week4');
        set_graphics(h,'xtick',1:16,'ytick',1:16, 'ax','square');
        axis square;
        subplot(234);
        colormap(hot); 
        imagesc(Mp12,scale); colorbar;
        title('Patients Week12');
        set_graphics(h,'xtick',1:16,'ytick',1:16, 'ax','square');
        axis square;
        subplot(235);
        colormap(hot); 
        imagesc(Mp24,scale); colorbar;
        title('Patients Week24');
        set_graphics(h,'xtick',1:16,'ytick',1:16, 'ax','square');
        subplot(236);
        colormap(hot); 
        imagesc(Mp52,scale); colorbar;
        title('Patients Week52');
        set_graphics(h,'xtick',1:16,'ytick',1:16, 'ax','square');
       
        
        s   = style_sheet(sty_2grp_ac,'leg',{'patients','control'},'leglocation','northoutside');
    
        x1  = nanmean(Dc.mat,1);
        
        DpW1 = getrow (Dp,Dp.week==1);
        x2 = nanmean(DpW1.mat,1);
        
        DpW4 = getrow (Dp,Dp.week==4);
        x4 = nanmean(DpW4.mat,1);
        
        DpW12 = getrow (Dp,Dp.week==12);
        x12 = nanmean(DpW12.mat,1);
        
        DpW24 = getrow (Dp,Dp.week==24);
        x24 = nanmean(DpW24.mat,1);
        
        DpW52 = getrow (Dp,Dp.week==52);
        x52 = nanmean(DpW52.mat,1);        
        
        h=make_figure;
        subplot(231);
        title('Week1');
        ylabel('correlation index')
        traceplot(1:25,[x1; x2],'split',[1 2]');
        subplot(232);
        traceplot(1:25,[x1; x4],'split',[1 2]');
        title('Week4');
        ylabel('correlation index')
        subplot(233);
        traceplot(1:25,[x1; x12],'split',[1 2]');
        title('Week12');
        ylabel('correlation index')
        subplot(234);
        traceplot(1:25,[x1; x24],'split',[1 2]');
        title('Week24');
        ylabel('correlation index')
        subplot(235);
        traceplot(1:25,[x1; x52],'split',[1 2]');
        title('Week52');
        ylabel('correlation index')
        xlabel('roi pairs')
        
        h=make_figure;
        subplot(231);
        scatterplot(x1',x2','regression','linear','printcorr')
        title('Week1');
        ylabel('control')
        xlabel('patients')
        axis square;
        subplot(232);
        scatterplot(x1',x4','regression','linear','printcorr')
        title('Week4');
        ylabel('control')
        xlabel('patients')
        axis square;
        subplot(233);
        scatterplot(x1',x12','regression','linear','printcorr')
        title('Week12');
        ylabel('control')
        xlabel('patients')
        axis square;
        subplot(234);
        scatterplot(x1',x24','regression','linear','printcorr')
        title('Week24');
        ylabel('control')
        xlabel('patients')
        axis square;
        subplot(235);
        scatterplot(x1',x52','regression','linear','printcorr')
        title('Week52');
        ylabel('control')
        xlabel('patients')
        axis square;
        
        h=make_figure;
        z1   = x1'-x2';
        z4   = x1'-x4';
        z12  = x1'-x12';
        z24  = x1'-x24';
        z52  = x1'-x52';
        subplot(231)
        plot(z1)
        title('Week1');
        ylabel('residuals')
        subplot(232)
        plot(z4)
        title('Week4');
        ylabel('residuals')
        subplot(233)
        plot(z12)
        title('Week12');
        ylabel('residuals')
        subplot(234)
        plot(z24)
        title('Week24');
        ylabel('residuals')
        subplot(235)
        plot(z52)
        title('Week52');
        ylabel('residuals')
        xlabel('pairs')
        
        keyboard;

        
    case 'plot_networkCorr_lesioned' % plot controls and patients (lesioned hemisphere only)
        D   = load(fullfile(rootDir,'rs_networkCorr_lesionedIntra.mat'));        
        
        %take only subgroup
        D   = getrow(D,(D.lesionType~=2)); %no cortical
        D   = getrow(D,(D.lesionType~=4)); %no mixed

            
%     only take chosen ones
%      subj = {'CU_2925','JHU_2531','UZ_2654','UZ_3239','UZ_3240','UZ_3241',...
%         'UZ_3243','UZ_3246','UZ_3247','UZ_3248','CU_1002','JHP_1001',...
%         'JHP_1002','JHP_1004','UZP_1001','UZP_1002','UZP_1004','UZP_1005',...
%         'UZP_1006','UZP_1007','UZP_1008'};
%      D = getrow(D,ismember(D.subj_name,subj));
%     
%     only exclude the bad ones 
% subj = {'CU_2697','JHU_2374','UZ_3228'};
%    D = getrow(D,~ismember(D.subj_name,subj));
        
        
        Dc  = getrow(D,D.control==1);
        Dp  = getrow(D,D.control~=1);
        
        % get avg matrices for each group
        Mc  = nanmean(Dc.M,1);
        Mc  = nonsym_squareform(Mc);
        
        
        DpW1 = getrow (Dp,Dp.week==1);
        Mp1 = nanmean(DpW1.M,1);
        Mp1 = nonsym_squareform(Mp1);
        DpW4 = getrow (Dp,Dp.week==4);
        Mp4 = nanmean(DpW4.M,1);
        Mp4 = nonsym_squareform(Mp4);
        DpW12 = getrow (Dp,Dp.week==12);
        Mp12 = nanmean(DpW12.M,1);
        Mp12 = nonsym_squareform(Mp12);
        DpW24 = getrow (Dp,Dp.week==24);
        Mp24 = nanmean(DpW24.M,1);
        Mp24 = nonsym_squareform(Mp24);
        DpW52 = getrow (Dp,Dp.week==52);
        Mp52 = nanmean(DpW52.M,1);
        Mp52 = nonsym_squareform(Mp52);
        
%         Mp = [];
%         for w=unique(Dp.week)'
%             Dw  = getrow(Dp,Dp.week==w);
%             Mpw  = nanmean(Dw.M,1);
%             Mpw  = nonsym_squareform(Mpw);
%             Mp = addstruct(Mp,Mpw);
%         end
        scale = [0 0.8];
        h=make_figure;
        subplot(231);
        colormap(hot); 
        imagesc(Mc,scale); colorbar;
        title('Controls');
        set_graphics(h,'xtick',1:16,'ytick',1:16, 'ax','square')
        subplot(232);
        colormap(hot); 
        imagesc(Mp1,scale); colorbar;
        title('Patients Week1');
        set_graphics(h,'xtick',1:16,'ytick',1:16, 'ax','square')
        subplot(233);
        colormap(hot); 
        imagesc(Mp4,scale); colorbar;
        title('Patients Week4');
        set_graphics(h,'xtick',1:16,'ytick',1:16, 'ax','square');
        axis square;
        subplot(234);
        colormap(hot); 
        imagesc(Mp12,scale); colorbar;
        title('Patients Week12');
        set_graphics(h,'xtick',1:16,'ytick',1:16, 'ax','square');
        axis square;
        subplot(235);
        colormap(hot); 
        imagesc(Mp24,scale); colorbar;
        title('Patients Week24');
        set_graphics(h,'xtick',1:16,'ytick',1:16, 'ax','square');
        subplot(236);
        colormap(hot); 
        imagesc(Mp52,scale); colorbar;
        title('Patients Week52');
        set_graphics(h,'xtick',1:16,'ytick',1:16, 'ax','square');
       
        
        s   = style_sheet(sty_2grp_ac,'leg',{'patients','control'},'leglocation','northoutside');
    
        h=make_figure;
        subplot(231);
        x1 = squareform(Mc);
        x2 = squareform(Mp1);
        title('Week1');
        ylabel('correlation index')
        traceplot(1:10,[x1; x2],'split',[1 2]');
        subplot(232);
        x4 = squareform(Mp4);
        traceplot(1:10,[x1; x4],'split',[1 2]');
        title('Week4');
        ylabel('correlation index')
        subplot(233);
        x12 = squareform(Mp12);
        traceplot(1:10,[x1; x12],'split',[1 2]');
        title('Week12');
        ylabel('correlation index')
        subplot(234);
        x24 = squareform(Mp24);
        traceplot(1:10,[x1; x24],'split',[1 2]');
        title('Week24');
        ylabel('correlation index')
        subplot(235);
        x52= squareform(Mp52);
        traceplot(1:10,[x1; x52],'split',[1 2]');
        title('Week52');
        ylabel('correlation index')
        xlabel('roi pairs')
        
        h=make_figure;
        subplot(231);
        scatterplot(x1',x2','regression','linear','printcorr')
        title('Week1');
        ylabel('control')
        xlabel('patients')
        axis square;
        subplot(232);
        x4 = squareform(Mp4);
        scatterplot(x1',x4','regression','linear','printcorr')
        title('Week4');
        ylabel('control')
        xlabel('patients')
        axis square;
        subplot(233);
        x12 = squareform(Mp12);
        scatterplot(x1',x12','regression','linear','printcorr')
        title('Week12');
        ylabel('control')
        xlabel('patients')
        axis square;
        subplot(234);
        x24 = squareform(Mp24);
        scatterplot(x1',x24','regression','linear','printcorr')
        title('Week24');
        ylabel('control')
        xlabel('patients')
        axis square;
        subplot(235);
        x52= squareform(Mp52);
        scatterplot(x1',x52','regression','linear','printcorr')
        title('Week52');
        ylabel('control')
        xlabel('patients')
        axis square;
        
        h=make_figure;
        z1   = x1'-x2';
        z4   = x1'-x4';
        z12  = x1'-x12';
        z24  = x1'-x24';
        z52  = x1'-x52';
        subplot(231)
        plot(z1)
        title('Week1');
        ylabel('residuals')
        subplot(232)
        plot(z4)
        title('Week4');
        ylabel('residuals')
        subplot(233)
        plot(z12)
        title('Week12');
        ylabel('residuals')
        subplot(234)
        plot(z24)
        title('Week24');
        ylabel('residuals')
        subplot(235)
        plot(z52)
        title('Week52');
        ylabel('residuals')
        xlabel('pairs')
        
        keyboard;        
    case 'plot_networkCorr_nonlesioned' % plot controls and patients (non lesioned hemisphere only)
        D   = load(fullfile(rootDir,'rs_networkCorr_nonlesionedIntra.mat'));        
        
        %take only subgroup
        D   = getrow(D,(D.lesionType~=2)); %no cortical
        D   = getrow(D,(D.lesionType~=4)); %no mixed

        Dc  = getrow(D,D.control==1);
        Dp  = getrow(D,D.control~=1);
        
        % get avg matrices for each group
        Mc  = nanmean(Dc.M,1);
        Mc  = nonsym_squareform(Mc);
        
        
        DpW1 = getrow (Dp,Dp.week==1);
        Mp1 = nanmean(DpW1.M,1);
        Mp1 = nonsym_squareform(Mp1);
        DpW4 = getrow (Dp,Dp.week==4);
        Mp4 = nanmean(DpW4.M,1);
        Mp4 = nonsym_squareform(Mp4);
        DpW12 = getrow (Dp,Dp.week==12);
        Mp12 = nanmean(DpW12.M,1);
        Mp12 = nonsym_squareform(Mp12);
        DpW24 = getrow (Dp,Dp.week==24);
        Mp24 = nanmean(DpW24.M,1);
        Mp24 = nonsym_squareform(Mp24);
        DpW52 = getrow (Dp,Dp.week==52);
        Mp52 = nanmean(DpW52.M,1);
        Mp52 = nonsym_squareform(Mp52);
        
%         Mp = [];
%         for w=unique(Dp.week)'
%             Dw  = getrow(Dp,Dp.week==w);
%             Mpw  = nanmean(Dw.M,1);
%             Mpw  = nonsym_squareform(Mpw);
%             Mp = addstruct(Mp,Mpw);
%         end
        scale = [0 0.8];
        h=make_figure;
        subplot(231);
        colormap(hot); 
        imagesc(Mc,scale); colorbar;
        title('Controls');
        set_graphics(h,'xtick',1:16,'ytick',1:16, 'ax','square')
        subplot(232);
        colormap(hot); 
        imagesc(Mp1,scale); colorbar;
        title('Patients Week1');
        set_graphics(h,'xtick',1:16,'ytick',1:16, 'ax','square')
        subplot(233);
        colormap(hot); 
        imagesc(Mp4,scale); colorbar;
        title('Patients Week4');
        set_graphics(h,'xtick',1:16,'ytick',1:16, 'ax','square');
        axis square;
        subplot(234);
        colormap(hot); 
        imagesc(Mp12,scale); colorbar;
        title('Patients Week12');
        set_graphics(h,'xtick',1:16,'ytick',1:16, 'ax','square');
        axis square;
        subplot(235);
        colormap(hot); 
        imagesc(Mp24,scale); colorbar;
        title('Patients Week24');
        set_graphics(h,'xtick',1:16,'ytick',1:16, 'ax','square');
        subplot(236);
        colormap(hot); 
        imagesc(Mp52,scale); colorbar;
        title('Patients Week52');
        set_graphics(h,'xtick',1:16,'ytick',1:16, 'ax','square');
       
        
        s   = style_sheet(sty_2grp_ac,'leg',{'patients','control'},'leglocation','northoutside');
    
        h=make_figure;
        subplot(231);
        x1 = squareform(Mc);
        x2 = squareform(Mp1);
        title('Week1');
        ylabel('correlation index')
        traceplot(1:10,[x1; x2],'split',[1 2]');
        subplot(232);
        x4 = squareform(Mp4);
        traceplot(1:10,[x1; x4],'split',[1 2]');
        title('Week4');
        ylabel('correlation index')
        subplot(233);
        x12 = squareform(Mp12);
        traceplot(1:10,[x1; x12],'split',[1 2]');
        title('Week12');
        ylabel('correlation index')
        subplot(234);
        x24 = squareform(Mp24);
        traceplot(1:10,[x1; x24],'split',[1 2]');
        title('Week24');
        ylabel('correlation index')
        subplot(235);
        x52= squareform(Mp52);
        traceplot(1:10,[x1; x52],'split',[1 2]');
        title('Week52');
        ylabel('correlation index')
        xlabel('roi pairs')
        
        h=make_figure;
        subplot(231);
        scatterplot(x1',x2','regression','linear','printcorr')
        title('Week1');
        ylabel('control')
        xlabel('patients')
        axis square;
        subplot(232);
        x4 = squareform(Mp4);
        scatterplot(x1',x4','regression','linear','printcorr')
        title('Week4');
        ylabel('control')
        xlabel('patients')
        axis square;
        subplot(233);
        x12 = squareform(Mp12);
        scatterplot(x1',x12','regression','linear','printcorr')
        title('Week12');
        ylabel('control')
        xlabel('patients')
        axis square;
        subplot(234);
        x24 = squareform(Mp24);
        scatterplot(x1',x24','regression','linear','printcorr')
        title('Week24');
        ylabel('control')
        xlabel('patients')
        axis square;
        subplot(235);
        x52= squareform(Mp52);
        scatterplot(x1',x52','regression','linear','printcorr')
        title('Week52');
        ylabel('control')
        xlabel('patients')
        axis square;
        
        h=make_figure;
        z1   = x1'-x2';
        z4   = x1'-x4';
        z12  = x1'-x12';
        z24  = x1'-x24';
        z52  = x1'-x52';
        subplot(231)
        plot(z1)
        title('Week1');
        ylabel('residuals')
        subplot(232)
        plot(z4)
        title('Week4');
        ylabel('residuals')
        subplot(233)
        plot(z12)
        title('Week12');
        ylabel('residuals')
        subplot(234)
        plot(z24)
        title('Week24');
        ylabel('residuals')
        subplot(235)
        plot(z52)
        title('Week52');
        ylabel('residuals')
        xlabel('pairs')
        
        keyboard;                
        
    case 'split_mvc'
        %D   = load(fullfile(rootDir,'rs_networkCorr_lesionedIntra.mat'));                
        %D   = load(fullfile(rootDir,'rs_networkCorr_nonlesionedIntra.mat'));                
        %D   = load(fullfile(rootDir,'rs_networkCorr.mat'));                
        
        
        for i=1:length(D.subj_name)
            m           = nonsym_squareform(D.M(i,:));
            D.mat(i,:)  = squareform(m);               %storing 10 pairs only
            %m           = m(:,1:5);
            %m           = m(6:10,:);
            D.mat(i,:)  = nonsym_squareform(m);
        end;
        
        %take only subgroup
        D   = getrow(D,(D.lesionType~=2)); %no cortical
        D   = getrow(D,(D.lesionType~=4)); %no mixed
        %sn  = {'UZ_3030','UZ_3226','UZ_3227'};
        %D   = getrow(D,~ismember(D.subj_name,sn));
        %D   = getrow(D,~(ismember(D.subj_name,'JHU_2713') & D.week==12));

        % 0. Avg controls (all) and patients (per week and according to severity)
        C           = tapply(D,{'subj_name'},{'mat','mean(x,1)','subset',D.control==1});

        D.mvcNorm   = -D.mvcNorm;                       % for some reason mvcNorm is negative (to be fixed)
        [x,sn]      = pivottable(D.subj_name,D.week,D.mvcNorm,'mean(x)','subset',D.control==0);
        x           = nanmean(x(:,1:2),2);
        
        weeks = [1 4 12 24 52];
        
        for w = 1:length(weeks)
            
        Dp          = getrow(D,D.control==0 & D.week==weeks(w) & ismember(D.subj_name,sn(~isnan(x))));
        Dp.sev      = zeros(length(Dp.subj_name),1);
        Dp.sev(ismember(Dp.subj_name,sn(x<0.3)))            = 1;    % mvc <0.3 at weeks 1-4
        Dp.sev(ismember(Dp.subj_name,sn(x>=0.3)))           = 1;    % mvc >=0.3 at weeks 1-4
        
        P           = tapply(Dp,{'sev'},{'mat','nanmean(x,1)'},...
                                        {'mvcNorm','nanmean(x)'});

        % 1. Compute correlations of
        %       - each control ens with left-out control ens mean
        %       - each control mm with left out control mm
        S = [];
        for i=1:length(C.subj_name)
            idxSubj     = i;
            idxLeftOut  = ~ismember(1:length(C.subj_name),i);

            rMat        = corr(C.mat(idxSubj,:)',mean(C.mat(idxLeftOut,:),1)');
            
            Si.subj_name    = C.subj_name(i);
            Si.rMat         = rMat;
            Si.mat          = C.mat(idxSubj,:);
            S           = addstruct(S,Si);
        end;

        % 2. Compute correlations of
        %       - each control ens with patient ens mean separately for each week
        %       - each control mm with patient mm mean separately for each week
        T = [];
        for i=1:length(C.subj_name)
            idxSubj     = i;

            rMat        = corr(P.mat',C.mat(idxSubj,:)');
            
            Ti.subj_name    = repmat(C.subj_name(i),length(rMat),1);
            Ti.sev          = P.sev;
            Ti.rMat         = rMat;
            T           = addstruct(T,Ti);
        end;
        varargout = {S,T};
        
        s   = style_sheet(sty_2grp_ac);
        subplot (2,3,w)
        myboxplot(T.sev,T.rMat,'split',T.sev,'CAT',s(1).CAT,s(1).PLOT{:});
        drawline(mean(S.rMat),'dir','horz','color','b');
        drawline(mean(S.rMat)+stderr(S.rMat),'dir','horz','color','b','linestyle','--');
        drawline(mean(S.rMat)-stderr(S.rMat),'dir','horz','color','b','linestyle','--');
        title(weeks(w))
        ylabel('pattern correlation')
        xlabel('severity') 
        
        fprintf('Weeks %d\n',weeks(w));
        fprintf('Sev 1\n');
        ttest(S.rMat,T.rMat(T.sev==1),2,'paired');
        fprintf('Sev 2\n');
        ttest(S.rMat,T.rMat(T.sev==2),2,'paired');
        fprintf('\n\n');
        end
    case 'split_mvc_hemis'
        %before you run it, check: median or mean? MVC or FM
        
        
         D   = load(fullfile(rootDir,'rs_networkCorr_lesionedIntra.mat'));                
         %D   = load(fullfile(rootDir,'rs_networkCorr_nonlesionedIntra.mat'));                               
        
        
        for i=1:length(D.subj_name)
            m           = nonsym_squareform(D.M(i,:));
            D.mat(i,:)  = squareform(m);               %storing 10 pairs only
        end;
        
        %take only subgroup
        D   = getrow(D,(D.lesionType~=2)); %no cortical
        D   = getrow(D,(D.lesionType~=4)); %no mixed

        % 0. Avg controls (all) and patients (per week and according to severity)
        %C           = tapply(D,{'subj_name'},{'mat','nanmean(x,1)','subset',D.control==1});
        C           = tapply(D,{'subj_name'},{'mat','nanmedian(x,1)','subset',D.control==1});
        
        %D.mvcNorm   = -D.mvcNorm;                       % for some reason mvcNorm is negative (to be fixed)
        %[x,sn]      = pivottable(D.subj_name,D.week,D.mvcNorm,'mean(x)','subset',D.control==0);
        %[x,sn]      = pivottable(D.subj_name,D.week,D.mvcNorm,'nanmedian(x)','subset',D.control==0);
        [x,sn]      = pivottable(D.subj_name,D.week,D.FM,'nanmedian(x)','subset',D.control==0);
        
        x           = nanmean(x(:,1:2),2);
        
        weeks = [1 4 12 24 52];
        
        for w = 1:length(weeks)
            
        Dp          = getrow(D,D.control==0 & D.week==weeks(w) & ismember(D.subj_name,sn(~isnan(x))));
        Dp.sev      = zeros(length(Dp.subj_name),1);
        
        %Dp.sev(ismember(Dp.subj_name,sn(x<0.3)))            = 1;    % mvc <0.3 at weeks 1-4
        %Dp.sev(ismember(Dp.subj_name,sn(x>=0.3)))           = 2;    % mvc >=0.3 at weeks 1-4
        
        Dp.sev(ismember(Dp.subj_name,sn(x<20)))            = 1;    % mvc <0.3 at weeks 1-4
        Dp.sev(ismember(Dp.subj_name,sn(x>=20)))           = 2;    % mvc >=0.3 at weeks 1-4
        
        %P           = tapply(Dp,{'sev'},{'mat','nanmean(x,1)'},...
        P           = tapply(Dp,{'sev'},{'mat','nanmedian(x,1)'},...
                                   {'FM','nanmedian(x)'});
                               %{'mvcNorm','nanmedian(x)'});
                                    
        % 1. Compute correlations of
        %       - each control ens with left-out control ens mean
        %       - each control mm with left out control mm
        S = [];
        for i=1:length(C.subj_name)
            idxSubj     = i;
            idxLeftOut  = ~ismember(1:length(C.subj_name),i);

            rMat        = corr(C.mat(idxSubj,:)',nanmedian(C.mat(idxLeftOut,:),1)');
            
            Si.subj_name    = C.subj_name(i);
            Si.rMat         = rMat;
            Si.mat          = C.mat(idxSubj,:);
            S           = addstruct(S,Si);
        end;

        % 2. Compute correlations of
        %       - each control ens with patient ens mean separately for each week
        %       - each control mm with patient mm mean separately for each week
        T = [];
        for i=1:length(C.subj_name)
            idxSubj     = i;

            rMat        = corr(P.mat',C.mat(idxSubj,:)');
            
            Ti.subj_name    = repmat(C.subj_name(i),length(rMat),1);
            Ti.sev          = P.sev;
            Ti.rMat         = rMat;
            T           = addstruct(T,Ti);
        end;
        varargout = {S,T};
        
        
        s   = style_sheet(sty_2grp_ac);
        subplot (2,3,w)
        myboxplot(T.sev,T.rMat,'split',T.sev,'CAT',s(1).CAT,s(1).PLOT{:});
        drawline(mean(S.rMat),'dir','horz','color','b');
        drawline(mean(S.rMat)+stderr(S.rMat),'dir','horz','color','b','linestyle','--');
        drawline(mean(S.rMat)-stderr(S.rMat),'dir','horz','color','b','linestyle','--');
        title(weeks(w))
        ylabel('pattern correlation')
        xlabel('severity') 
        
        fprintf('Weeks %d\n',weeks(w));
        fprintf('Sev 1\n');
        ttest(S.rMat,T.rMat(T.sev==1),2,'paired');
        fprintf('Sev 2\n');
        ttest(S.rMat,T.rMat(T.sev==2),2,'paired');
        fprintf('\n\n');
        end
    case 'split_mvc_inter'
    %first check if median or mean and FM or mvc and split or no split
        
    D   = load(fullfile(rootDir,'rs_networkCorr.mat')); 
    
    for i=1:length(D.subj_name)
            m           = nonsym_squareform(D.M(i,:));
            m           = m(:,1:5);
            m           = m(6:10,:);
            D.mat(i,:)  = nonsym_squareform(m);
    end;

    %take only subgroup
        D   = getrow(D,(D.lesionType~=2)); %no cortical
        D   = getrow(D,(D.lesionType~=4)); %no mixed

        % 0. Avg controls (all) and patients (per week and according to severity)
        C           = tapply(D,{'subj_name'},{'mat','nanmedian(x,1)','subset',D.control==1});

        %D.mvcNorm   = -D.mvcNorm;                       % for some reason mvcNorm is negative (to be fixed)
        %[x,sn]      = pivottable(D.subj_name,D.week,D.mvcNorm,'nanmedian(x)','subset',D.control==0);
        [x,sn]      = pivottable(D.subj_name,D.week,D.FM,'nanmedian(x)','subset',D.control==0);
        
        x           = nanmean(x(:,1:2),2);
        
        weeks = [1 4 12 24 52];
        
        for w = 1:length(weeks)
            
        Dp          = getrow(D,D.control==0 & D.week==weeks(w) & ismember(D.subj_name,sn(~isnan(x))));
        Dp.sev      = zeros(length(Dp.subj_name),1);
        Dp.sev(ismember(Dp.subj_name,sn(x<20)))            = 1;    % mvc <0.3 at weeks 1-4
       %change accordingly
        Dp.sev(ismember(Dp.subj_name,sn(x>=20)))           = 2;    % mvc >=0.3 at weeks 1-4
        
        P           = tapply(Dp,{'sev'},{'mat','nanmedian(x,1)'},...
                                        {'mvcNorm','nanmedian(x)'});

        % 1. Compute correlations of
        %       - each control ens with left-out control ens mean
        %       - each control mm with left out control mm
        S = [];
        for i=1:length(C.subj_name)
            idxSubj     = i;
            idxLeftOut  = ~ismember(1:length(C.subj_name),i);

            rMat        = corr(C.mat(idxSubj,:)',nanmedian(C.mat(idxLeftOut,:),1)');
            
            Si.subj_name    = C.subj_name(i);
            Si.rMat         = rMat;
            Si.mat          = C.mat(idxSubj,:);
            S           = addstruct(S,Si);
        end;

        % 2. Compute correlations of
        %       - each control ens with patient ens mean separately for each week
        %       - each control mm with patient mm mean separately for each week
        T = [];
        for i=1:length(C.subj_name)
            idxSubj     = i;

            rMat        = corr(P.mat',C.mat(idxSubj,:)');
            
            Ti.subj_name    = repmat(C.subj_name(i),length(rMat),1);
            Ti.sev          = P.sev;
            Ti.rMat         = rMat;
            T           = addstruct(T,Ti);
        end;
        varargout = {S,T};
        
        
        s   = style_sheet(sty_2grp_ac);
        subplot (2,3,w)
        myboxplot(T.sev,T.rMat,'split',T.sev,'CAT',s(1).CAT,s(1).PLOT{:});   %ask Naveed
        %myboxplot(T.sev,T.rMat,'split',T.sev,'CAT',s(1).CAT,s(1).PLOT{:});
        drawline(mean(S.rMat),'dir','horz','color','b');
        drawline(mean(S.rMat)+stderr(S.rMat),'dir','horz','color','b','linestyle','--');
        drawline(mean(S.rMat)-stderr(S.rMat),'dir','horz','color','b','linestyle','--');
        title(weeks(w))
        ylabel('pattern correlation')
        %xlabel('severity') 
        
        fprintf('Weeks %d\n',weeks(w));
        fprintf('Sev 1\n');
        ttest(S.rMat,T.rMat(T.sev==1),2,'paired');
         fprintf('Sev 2\n');
         ttest(S.rMat,T.rMat(T.sev==2),2,'paired');
         fprintf('\n\n');
        end 
    keyboard;   
    case 'split_mvcMagnitude'
        D   = load(fullfile(rootDir,'rs_networkCorr_lesionedIntra.mat'));                
        % D   = load(fullfile(rootDir,'rs_networkCorr_nonlesionedIntra.mat'));                
%         D   = load(fullfile(rootDir,'rs_networkCorr.mat'));                
        
        for i=1:length(D.subj_name)
            m           = nonsym_squareform(D.M(i,:));
            D.rMag(i,1) = mean(squareform(m));
        end;
        
        %take only subgroup
        D   = getrow(D,(D.lesionType~=2)); %no cortical
        D   = getrow(D,(D.lesionType~=4)); %no mixed

        % 0. Avg controls (all) and patients (per week and according to severity)
        C           = tapply(D,{'subj_name'},{'rMag','mean(x)','subset',D.control==1});

        D.mvcNorm   = -D.mvcNorm;                       % for some reason mvcNorm is negative (to be fixed)
        [x,sn]      = pivottable(D.subj_name,D.week,D.mvcNorm,'mean(x)','subset',D.control==0);
        x           = nanmean(x(:,1:2),2);
        
        weeks = [1 4 12 24 52];
        
        for w = 1:length(weeks)
            
        Dp          = getrow(D,D.control==0 & D.week==weeks(w) & ismember(D.subj_name,sn(~isnan(x))));
        Dp.sev      = zeros(length(Dp.subj_name),1);
        Dp.sev(ismember(Dp.subj_name,sn(x<0.3)))            = 1;    % mvc <0.3 at weeks 1-4
        Dp.sev(ismember(Dp.subj_name,sn(x>=0.3)))           = 2;    % mvc >=0.3 at weeks 1-4
        
        P           = tapply(Dp,{'subj_name','sev'},{'rMag','nanmean(x,1)'},...
                                        {'mvcNorm','nanmean(x)'});

        s   = style_sheet(sty_2grp_ac);
        subplot (2,3,w)
        myboxplot(P.sev,P.rMag,'split',P.sev,'CAT',s(1).CAT,s(1).PLOT{:});
        drawline(mean(C.rMag),'dir','horz','color','b');
        drawline(mean(C.rMag)+stderr(C.rMag),'dir','horz','color','b','linestyle','--');
        drawline(mean(C.rMag)-stderr(C.rMag),'dir','horz','color','b','linestyle','--');
        title(weeks(w))
        ylabel('pattern correlation')
        xlabel('severity') 
        
        fprintf('Weeks %d\n',weeks(w));
        fprintf('Controls vs Sev 1\n');
        ttest(C.rMag,P.rMag(P.sev==1),2,'independent');
        fprintf('Control vs Sev 2\n');
        ttest(C.rMag,P.rMag(P.sev==2),2,'independent');
        fprintf('Sev 1 vs Sev 2\n');
        ttest(P.rMag(P.sev==1),P.rMag(P.sev==2),2,'independent');
        fprintf('\n\n');
        end
        
    case 'plot_networkCorr_nonlesioned_cortical' % plot controls and patients (non lesioned hemisphere only)
        D   = load(fullfile(rootDir,'rs_networkCorr_nonlesionedIntra.mat'));        
        
        %take only subgroup
        D   = getrow(D,(D.lesionType~=3)); %no subcortical
        sn  = {'UZ_3030','UZ_3226','UZ_3227'};
        D   = getrow(D,~ismember(D.subj_name,sn));
        D   = getrow(D,~(ismember(D.subj_name,'JHU_2713') & D.week==12));

        Dc  = getrow(D,D.control==1);
        Dp  = getrow(D,D.control~=1);
        
        % get avg matrices for each group
        Mc  = nanmean(Dc.M,1);
        Mc  = nonsym_squareform(Mc);
        
        
        DpW1 = getrow (Dp,Dp.week==1);
        Mp1 = nanmean(DpW1.M,1);
        Mp1 = nonsym_squareform(Mp1);
        DpW4 = getrow (Dp,Dp.week==4);
        Mp4 = nanmean(DpW4.M,1);
        Mp4 = nonsym_squareform(Mp4);
        DpW12 = getrow (Dp,Dp.week==12);
        Mp12 = nanmean(DpW12.M,1);
        Mp12 = nonsym_squareform(Mp12);
        DpW24 = getrow (Dp,Dp.week==24);
        Mp24 = nanmean(DpW24.M,1);
        Mp24 = nonsym_squareform(Mp24);
        DpW52 = getrow (Dp,Dp.week==52);
        Mp52 = nanmean(DpW52.M,1);
        Mp52 = nonsym_squareform(Mp52);
        
%         Mp = [];
%         for w=unique(Dp.week)'
%             Dw  = getrow(Dp,Dp.week==w);
%             Mpw  = nanmean(Dw.M,1);
%             Mpw  = nonsym_squareform(Mpw);
%             Mp = addstruct(Mp,Mpw);
%         end
        scale = [0 0.8];
        h=make_figure;
        subplot(231);
        colormap(hot); 
        imagesc(Mc,scale); colorbar;
        title('Controls');
        set_graphics(h,'xtick',1:16,'ytick',1:16, 'ax','square')
        subplot(232);
        colormap(hot); 
        imagesc(Mp1,scale); colorbar;
        title('Patients Week1');
        set_graphics(h,'xtick',1:16,'ytick',1:16, 'ax','square')
        subplot(233);
        colormap(hot); 
        imagesc(Mp4,scale); colorbar;
        title('Patients Week4');
        set_graphics(h,'xtick',1:16,'ytick',1:16, 'ax','square');
        axis square;
        subplot(234);
        colormap(hot); 
        imagesc(Mp12,scale); colorbar;
        title('Patients Week12');
        set_graphics(h,'xtick',1:16,'ytick',1:16, 'ax','square');
        axis square;
        subplot(235);
        colormap(hot); 
        imagesc(Mp24,scale); colorbar;
        title('Patients Week24');
        set_graphics(h,'xtick',1:16,'ytick',1:16, 'ax','square');
        subplot(236);
        colormap(hot); 
        imagesc(Mp52,scale); colorbar;
        title('Patients Week52');
        set_graphics(h,'xtick',1:16,'ytick',1:16, 'ax','square');
       
        
        s   = style_sheet(sty_2grp_ac,'leg',{'patients','control'},'leglocation','northoutside');
    
        h=make_figure;
        subplot(231);
        x1 = squareform(Mc);
        x2 = squareform(Mp1);
        title('Week1');
        ylabel('correlation index')
        traceplot(1:10,[x1; x2],'split',[1 2]');
        subplot(232);
        x4 = squareform(Mp4);
        traceplot(1:10,[x1; x4],'split',[1 2]');
        title('Week4');
        ylabel('correlation index')
        subplot(233);
        x12 = squareform(Mp12);
        traceplot(1:10,[x1; x12],'split',[1 2]');
        title('Week12');
        ylabel('correlation index')
        subplot(234);
        x24 = squareform(Mp24);
        traceplot(1:10,[x1; x24],'split',[1 2]');
        title('Week24');
        ylabel('correlation index')
        subplot(235);
        x52= squareform(Mp52);
        traceplot(1:10,[x1; x52],'split',[1 2]');
        title('Week52');
        ylabel('correlation index')
        xlabel('roi pairs')
        
        h=make_figure;
        subplot(231);
        scatterplot(x1',x2','regression','linear','printcorr')
        title('Week1');
        ylabel('control')
        xlabel('patients')
        axis square;
        subplot(232);
        x4 = squareform(Mp4);
        scatterplot(x1',x4','regression','linear','printcorr')
        title('Week4');
        ylabel('control')
        xlabel('patients')
        axis square;
        subplot(233);
        x12 = squareform(Mp12);
        scatterplot(x1',x12','regression','linear','printcorr')
        title('Week12');
        ylabel('control')
        xlabel('patients')
        axis square;
        subplot(234);
        x24 = squareform(Mp24);
        scatterplot(x1',x24','regression','linear','printcorr')
        title('Week24');
        ylabel('control')
        xlabel('patients')
        axis square;
        subplot(235);
        x52= squareform(Mp52);
        scatterplot(x1',x52','regression','linear','printcorr')
        title('Week52');
        ylabel('control')
        xlabel('patients')
        axis square;
        
        h=make_figure;
        z1   = x1'-x2';
        z4   = x1'-x4';
        z12  = x1'-x12';
        z24  = x1'-x24';
        z52  = x1'-x52';
        subplot(231)
        plot(z1)
        title('Week1');
        ylabel('residuals')
        subplot(232)
        plot(z4)
        title('Week4');
        ylabel('residuals')
        subplot(233)
        plot(z12)
        title('Week12');
        ylabel('residuals')
        subplot(234)
        plot(z24)
        title('Week24');
        ylabel('residuals')
        subplot(235)
        plot(z52)
        title('Week52');
        ylabel('residuals')
        xlabel('pairs')
        
        keyboard;            
        
    case 'otherwise'
        disp('No such case...');
        
%     case 'PP_additional'        % additional variables added to preprocessed data to make analysis easier
%         %D   = load(fullfile(ppDir,'rs_preprocess.mat'));
%         D       = load('rs_preprocess.mat');
%         
%         % 1. save fisher z transform r values
%         D.fz_r  = fisherz(D.r);
%         
%         % 2. add location for region side (left or right hem)
%         %   - 1 is left hem, 2 is right
%         D.hem           = regSide(D.rpair);
%         
%         % 3. add a variable call hl which is
%         %   - 0 if interhem, 1 if intra (left hem), 2 if intra (right hem)
%         D.hl            = mean(D.hem,2);
%         D.hl(D.hl==1.5) = 0;
%         
%         % 4. add a variable called regType which is
%         %   - 1 for S1, 2 for M1, 3 for PMd etc ...
%         D.regType   = regType(D.rpair);
%         
%         % 4. save structure
%         save('rs_preprocess.mat','-struct','D');
% %         keyboard;
%         
% %         % 5. How to use these variables
% %         % Example 1. Get m1-m1 correlations for controls
% %         T = getrow(D,D.control==1 & D.regType(:,1)==2 & D.regType(:,2)==2);
% %         % Example 2. Get intrahem m1-s1 and s1-m1 correlations for controls
% %         T = getrow(D,D.control==1 & D.regType(:,1)==1 & D.regType(:,2)==2);
% %         % Example 3. Get intrahem m1-s1 and s1-m1 correlations for controls
% %         % (only interhemispheric)
% %         T = getrow(D,D.control==1 & D.hl==0 & D.regType(:,1)==1 & D.regType(:,2)==2);
% %         % Example 4. Get intrahem m1-s1 and s1-m1 correlations for controls
% %         % (only intrahemispheric)
% %         T = getrow(D,D.control==1 & D.hl~=0 & D.regType(:,1)==1 & D.regType(:,2)==2);
% %         % Example 5. Get all lesioned intrahemispheric correlations
% %         T = getrow(D,D.control==0 & D.lesionSide==D.hl);    
% %         
% %         % Example 6. 
% %         T = getrow(D,D.control==0 & D.lesionSide==D.hl);    % Get lesioned hemisphere
% %         S = getrow(D,D.control==1 & D.hl~=0);                         % Get controls
% %         A = addstruct(T,S);
% % %         unique(A.hl);         % sanity check to see if hl is not 0 
% %         lineplot(A.week,A.r,'split',A.control==1,'subset',A.regType(:,1)==1 & A.regType(:,2)==2,...
% %                             'plotfcn','nanmean','leg','auto','style_thickline');

    otherwise
        disp('no such case');
end;