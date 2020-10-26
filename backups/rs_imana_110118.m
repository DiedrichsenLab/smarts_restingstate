 function varargout=rs_imana_110118(what,varargin)
% MB - analysis script for the SMARTS resting state data
% Open questions:
% - r or fisher z
% - dominant = nonlesion controls???


% rootDir             = '/Volumes/External/data/smarts';
rootDir             = '/Volumes/Naveed/data/smarts';
% rootDir             = '/Users/naveed/Documents/data/smarts';
% rootDir             = '/Users/naveed/Dropbox/Personal/N&M/rsfMRI_results';
% rootDir             = '/Users/meret/Dropbox/N&M/rsfMRI_results/resting_state';
%rootDir             = '/Volumes/Naveed/data/smarts';
behDir              = [rootDir '/bedside/analysis'];
rsDir               = [rootDir '/fmri/restingstate_imaging'];
roiDir              = [rootDir '/fmri/RegionOfInterest'];
% ppDir               = [rootDir '/new/'];
% ppDir               = [rsDir '/preprocess'];
analysisDir         = [rsDir '/preprocess'];
ppDir               = [rootDir];
preZipDir           = [rootDir '/fmri/restingstate_imaging/pre'];
statsDir            = [rootDir '/R_new/'];

roiPark2011         = [rsDir '/RegionOfInterest-Park2011'];
%rDir                = [rootDir, '/new/R'];

% color_scheme
colours             = num2cell(parula(20),2)';
sty_2grp_ac         = colours([8,12]);   % e.g. group 1 vs group 2
sty_2grp_cp         = colours([5,15]);
sty_multiple        = colours([1,3,6]);

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
        % 0. check which subjects have the raw and preprocessed files
        D = rs_preprocess_ana('check_data_exists');
        D = getrow(D,D.existPPImg & D.existRawImg);
        D = getrow(D,D.lesionside==0);      % get only controls

        % need to estimate rel. over the two preprocessing techniques
        %   pre:        simple pipeline
        %   preproc:    ICA-based noise removal pipeline
        type = {'pre_Bold_Rest','preproc_Bold_Rest'};          
        roi  = [2 10];  % left and right M1

        S = [];
        for i=1:length(D.SN)
            Di = getrow(D,i);       % get i-th session
            fprintf('%s - %s\n',Di.subj_name{1},Di.Week{1});
            
            % load ROI definitions for that subjection
            R   = load(fullfile(roiDir,sprintf('%s_regions_all.mat',Di.subj_name{1})));
            
            % loop over two preprocessing types
            for p=1:length(type)
                fdir    = fullfile(rsDir,Di.subj_name{1},Di.Week{1});
                fname   = sprintf('%s_%s_%s.nii',type{p},Di.subj_name{1},Di.Week{1});
                
                if exist(fullfile(fdir,fname),'file')
                    v   = spm_vol(fullfile(fdir,fname));        % load volumes from rs data
                    ts  = region_getdata(v,R.R(roi));           % get time series for M1-left and M1-right

                    % calculate corr for time series
                    ts1 = mean(ts{1},2);
                    ts2 = mean(ts{2},2);
                    r   = corr(ts1,ts2);
                else
                    r   = nan;
                end;
                
                % write out data to structure
                Si              = Di;
                Si.ptype        = p;
                Si.pname        = type(p);
                Si.r            = r;                % M1-M1 correlation
                S               = addstruct(S,Si);
            end;
        end;
        varargout = {S};
        save(fullfile(analysisDir,'rs_compare_preprocess.mat'),'-struct','S');
        
    case 'PP_rawTimeSeries'             % Extraction of fisher-z correlations from different regions of interest from the freesurfer atlas
        pre = {'preproc_Bold_Rest'};
        lesionType = {'none','cortical','subcortical','mixed'}; %1,2,3,4
        roi  = combnk(1:16,2); % all pairs
        
        % 0. check which subjects have the raw and preprocessed files
        D = rs_preprocess_ana('check_data_exists');
        D = getrow(D,D.existPPImg & D.existRawImg);
        
%         D   = dload(fullfile(rsDir,'patient_list.txt'));
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
                fDir= fullfile(rsDir,D.subj_name{i},D.Week{i});
                v   = spm_vol(fullfile(fDir,fname));       % load volumes from rs data
                ts  = region_getdata(v,R.R);                % get time series for all ROIs
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
                    Si.mvc          = Bj.mvc;
                    Si.FM           = Bj.FM;
                    Si.ARAT         = Bj.ARAT;
                    Si.age          = Bj.age;
                    Si.handedness   = Bj.handedness;
                    Si.lesion_PrCG  = SSi.lesion_PrCG(1);
                    Si.lesion_PrCG_W= SSi.lesion_PrCG_W(1);
                else
                    Si.ensOverall   = nan;
                    Si.mmOverall    = nan;
                    Si.mvcNorm      = nan;
                    Si.mvc          = nan;
                    Si.FM           = nan;
                    Si.ARAT         = nan;
                    Si.age          = nan;
                    Si.handedness   = nan;
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
        
        % 0. check which subjects have the raw and preprocessed files
        D = rs_preprocess_ana('check_data_exists');
        D = getrow(D,D.existPPImg & D.existRawImg);
        
%         D   = dload(fullfile(ppDir,'patient_list.txt'));
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
                fDir= fullfile(rsDir,D.subj_name{i},D.Week{i});
                v   = spm_vol(fullfile(fDir,fname));       % load volumes from rs data                
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
                    Si.mvc          = Bj.mvc;
                    Si.FM           = Bj.FM;
                    Si.ARAT         = Bj.ARAT;                    
                    Si.lesion_PrCG  = SSi.lesion_PrCG(1);
                    Si.lesion_PrCG_W= SSi.lesion_PrCG_W(1);
                else
                    Si.ensOverall   = nan;
                    Si.mmOverall    = nan;
                    Si.mvcNorm      = nan;
                    Si.mvc          = nan;
                    Si.ARAT         = nan;                    
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
    case 'PP_removePair'                % some participants have bad correlation pairs
        D = varargin{1};
        
        idx         = (D.r==0 | D.r==1 | D.r==-1);  % edge cases for bad correlation pairs
        D.r(idx)    = nan;
        D.fz_r      = fisherz(D.r);
        
        varargout = {D};        
    case 'PP_studySubset'               % get controls and patients (subcorticals only, more than 1 measurement session)
        D    = load(fullfile(ppDir,'rs_preprocess.mat'));
%         D    = load(fullfile(ppDir,'rs_preprocess_splithalf.mat'));
        
        % 1. get controls 
        Dc      = getrow(D,D.control==1);
        [x,sn]  = pivottable(Dc.subj_name,[],Dc.week,'length(unique(x))');
        Dc      = getrow(Dc,ismember(Dc.subj_name,sn(x>1)));
        
        % remove bad control
        idx     = ismember(Dc.subj_name,'CUP_1002') & Dc.week==24;
        Dc      = getrow(Dc,~idx);
        
        Dp      = getrow(D,D.lesionType==3);
        [x,sn]  = pivottable(Dp.subj_name,[],Dp.week,'length(unique(x))');
        Dp      = getrow(Dp,ismember(Dp.subj_name,sn(x>1)));
        
        T = addstruct(Dc,Dp);
        varargout = {T};
        
    case 'NC_getInterHemPattern'        % get pattern of correlations for pairs between hemispheres
        % 0. Load data
        D       = rs_imana_110118('PP_studySubset');
        D       = rs_imana_110118('PP_removePair',D);

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
                
                % flip hemispheres (left is non-lesioned, right is
                % lesioned)
                if mean(Dw.lesionSide)==1
                    M = M';     % flip to make correlation between pairs for non-lesioned first
                end;
                
                % add data to new structure
                Si      = getrow(Dw,1);
                Si      = rmfield(Si,{'rType','rpair','r'});
                Si.M    = nonsym_squareform(M);
                S       = addstruct(S,Si);
            end;
        end;
        save(fullfile(ppDir,'nc_interhem.mat'),'-struct','S'); 
    case 'NC_getIntraLesionedHemPattern'      % get pattern of correlations for pairs within lesioned hemisphere
        % 0. Load data
        D       = rs_imana_110118('PP_studySubset');
        D       = rs_imana_110118('PP_removePair',D);

        % 0. Obtain required data subset        
        %       - exclude unwanted ROIs
        inclROI = [1:5];        %no parietal or v1 rois         
        idx = ismember(D.regType(:,1),inclROI) & ismember(D.regType(:,2),inclROI);        
        D   = getrow(D,idx);        
        
        %       - get intrahemispheric fisherz correlations
        D   = getrow(D,D.hl~=0);
        Dc  = getrow(D,D.control==1);
        Dp  = getrow(D,D.control==0 & D.lesionSide==D.hl);
        clear D
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
        save(fullfile(ppDir,'nc_intrahem_lesioned.mat'),'-struct','S');
    case 'NC_getIntraNonlesionedHemPattern'   % get pattern of correlations for pairs within lesioned hemisphere
        % 0. Load data
        D       = rs_imana_110118('PP_studySubset');
        D       = rs_imana_110118('PP_removePair',D);

        % 0. Obtain required data subset        
        %       - exclude unwanted ROIs
        inclROI = [1:5];        %no parietal or v1 rois         
        idx = ismember(D.regType(:,1),inclROI) & ismember(D.regType(:,2),inclROI);        
        D   = getrow(D,idx);        
        
        %       - get intrahemispheric fisherz correlations
        D   = getrow(D,D.hl~=0);
        Dc  = getrow(D,D.control==1);
        Dp  = getrow(D,D.control==0 & D.lesionSide~=D.hl);
        clear D
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
        save(fullfile(ppDir,'nc_intrahem_nonlesioned.mat'),'-struct','S');         
    case 'NC_getThemALL'                      % get pattern of correlations for pairs inter/intra
              
       D3 = load(fullfile(ppDir,'nc_intrahem_lesioned.mat')); 
       D2 = load(fullfile(ppDir,'nc_intrahem_nonlesioned.mat')); 
       D1 = load(fullfile(ppDir,'nc_interhem.mat')); 
       
       S = D1;
       S.M =[D1.M, D2.M, D3.M];       
       save(fullfile(ppDir,'nc_getThemall.mat'),'-struct','S'); 
        
    case 'NC_getInterHemPatternSH'      % get pattern of correlations for pairs between hemispheres (split half)
        % 0. Load data
        D       = rs_imana_110118('PP_studySubset'); 
        D       = rs_imana_110118('PP_removePair',D);

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
                    
                    % flip hemispheres (left is non-lesioned, right is
                    % lesioned)
                    if mean(Dsp.lesionSide)==1
                        M(:,:,1) = M(:,:,1)';
                        M(:,:,2) = M(:,:,2)';
                        % flip to make correlation between pairs for non-lesioned first
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
        save(fullfile(ppDir,'nc_interhem_splithalf.mat'),'-struct','S'); 
    case 'NC_getIntraHemLesionedSH'     % get pattern of correlations for pairs between hemispheres (split half)
        % 0. Load data
        D       = rs_imana_110118('PP_studySubset'); 
        D       = rs_imana_110118('PP_removePair',D);

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
        save(fullfile(ppDir,'nc_intrahem_lesioned_splithalf.mat'),'-struct','S');   
    case 'NC_getIntraHemNonlesionedSH'  % get pattern of correlations for pairs between hemispheres (split half)
        % 0. Load data
        D       = rs_imana_110118('PP_studySubset'); 
        D       = rs_imana_110118('PP_removePair',D);

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
        save(fullfile(ppDir,'nc_intrahem_nonlesioned_splithalf.mat'),'-struct','S');           
    case 'NC_getThemALLSH'
     D3 = load(fullfile(ppDir,'nc_intrahem_lesioned_splithalf.mat')); 
     D2 = load(fullfile(ppDir,'nc_intrahem_nonlesioned_splithalf.mat')); 
     D1 = load(fullfile(ppDir,'nc_interhem_splithalf.mat')); 
       
       S = D1;
       S.M1=[D1.M1, D2.M1, D3.M1];  
       S.M2 =[D1.M2, D2.M2, D3.M2];  
       
        T   = [];
        % Loop over subj
        for sn=unique(S.subj_name)'
            Ds = getrow(S,strcmp(S.subj_name,sn));
            
            % Loop over week
            for w=unique(Ds.week)'
                Dw  = getrow(Ds,Ds.week==w);
                fprintf('%s - W%d\n',Dw.subj_name{1},Dw.week(1));
                
                % add data to new structure
                Si          = getrow(Dw,1);
                Si.M1       = Dw.M1;
                Si.M2       = Dw.M2;
                Si.sh_r     = corr(Dw.M1',Dw.M2');
                Si.sh_fz_r  = fisherz(Si.sh_r);
                T           = addstruct(T,Si);
            end;
        end;
       
       save(fullfile(ppDir,'nc_getThemall_splithalf.mat'),'-struct','T');    
        
    case 'REPLICATION_golestani2013_getRaw'
        pre = {'preproc_Bold_Rest'};
        lesionType = {'none','cortical','subcortical','mixed'}; %1,2,3,4
        roi  = combnk(1:16,2); % all pairs
        inclROI = [1 2 9 10];
        
        
        % 0. check which subjects have the raw and preprocessed files
        D = rs_preprocess_ana('check_data_exists');
        D = getrow(D,D.existPPImg & D.existRawImg);
        
        % remove bad control session
        idx = ismember(D.subj_name,'CUP_1002') & D.week==24;
        D   = getrow(D,~idx);
        
%         D   = dload(fullfile(ppDir,'patient_list.txt'));
        B   = load(fullfile(behDir,'ens_alldat_mvc.mat'));
        SS  = dload(fullfile(rootDir,'bedside','data_clean','subject_list.txt'));
        L   = dload(fullfile(rootDir,'DTI','lesion_analysis','lesion.dat'));
        
        SS.subj_name    = strcat(SS.Centre,'_',num2str(SS.ID));
        L.subj_name     = strcat(L.Centre,'_',num2str(L.ID));
        
        S = [];
        for i=1:length(D.SN)          
            if D.SN(i) == 186
                continue;
            end;
            
            fname   = strcat(pre{1},'_',D.Centre{i},'_',num2str(D.ID(i)),'_',D.Week{i},'.nii');
            sn      = strcat(D.Centre{i},'_',num2str(D.ID(i)));
            fprintf('%d.\t%s\n',D.SN(i),fname);
            
            regFile = fullfile(roiDir,sprintf('%s_regions_all.mat',sn));
            if ~exist(regFile)
                fprintf('No region definition for this subject...\n');
            else
                R   = load(fullfile(roiDir,sprintf('%s_regions_all.mat',sn))); % load subject region definition
                R.R = R.R(inclROI);
                fDir= fullfile(rsDir,D.subj_name{i},D.Week{i});
                v   = spm_vol(fullfile(fDir,fname));       % load volumes from rs data
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
                    Si.mvc          = Bj.mvc;
                    Si.FM           = Bj.FM;
                    Si.ARAT         = Bj.ARAT;
                    Si.lesion_PrCG  = SSi.lesion_PrCG(1);
                    Si.lesion_PrCG_W= SSi.lesion_PrCG_W(1);
                else
                    Si.ensOverall   = nan;
                    Si.mmOverall    = nan;
                    Si.mvcNorm      = nan;
                    Si.mvc          = nan;
                    Si.FM           = nan;
                    Si.ARAT         = nan;
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

                % estimate left and right combined M1/S1 regions
                ts1 = [ts{1},ts{2}];    % left hemisphere (M1/S1)
                ts2 = [ts{3},ts{4}];    % right hemisphere (M1/S1)
                
                % removing all voxels with no signal
                ts1 = ts1(:,all(ts1,1));    
                ts2 = ts2(:,all(ts2,1));    

                refTS   = []; % reference time-series (lesioned hem for patients, left hem for controls)
                compTS  = [];
                if (Si.control==1)
                    refTS   = ts1;
                    compTS  = ts2;
                else
                    if (Si.lesionSide==1)
                        refTS   = ts1;
                        compTS  = ts2;
                    else
                        refTS   = ts2;
                        compTS  = ts1;
                    end;
                end;
                
                % estimate pair-wise correlations
                %   - fisherz transform (distribution exists around -3 and 3)
                %   - remove inf
                r       = corr(refTS);
                idx     = any(~isnan(r),1);
                refFZ   = fisherz(1-squareform(1-r(idx,idx)));
%                 refFZ   = refFZ(refFZ>-3 & refFZ<3);
                
                r       = corr(refTS,compTS);
                compFZ  = fisherz(r(~isnan(r)));
                compFZ  = compFZ(:);
                compFZ  = compFZ(compFZ>-3 & compFZ<3);
                
                Si.compFZ   = mean(compFZ);
                Si.refFZ    = mean(refFZ);
%                 Si.refCon   = mean(Si.compFZ)./mean(Si.refFZ);
                
                S           = addstruct(S,Si);
            end;
        end;
        varargout = {S};
        save(fullfile(ppDir,'golestani2013_preprocess.mat'),'-struct','S');     
    case 'REPLICATION_park2011_getCorr'     
        pre = {'preproc_Bold_Rest'};
        lesionType = {'none','cortical','subcortical','mixed'}; %1,2,3,4
        roi  = combnk(1:16,2); % all pairs
        inclROI = [1 2 3 4];
        
        
        % 0. check which subjects have the raw and preprocessed files
        D = rs_preprocess_ana('check_data_exists');
        D = getrow(D,D.existPPImg & D.existRawImg);

        % remove bad control session
        idx = ismember(D.subj_name,'CUP_1002') & D.week==24;
        D   = getrow(D,~idx);        
        
%         D   = dload(fullfile(ppDir,'patient_list.txt'));
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
            
            regFile = fullfile(roiPark2011,sprintf('%s_regions_all.mat',sn));
            if ~exist(regFile)
                fprintf('No region definition for this subject...\n');
            else
                R   = load(fullfile(roiPark2011,sprintf('%s_regions_all.mat',sn))); % load subject region definition
                R.R = R.R(inclROI);
                fDir= fullfile(rsDir,D.subj_name{i},D.Week{i});
                v   = spm_vol(fullfile(fDir,fname));       % load volumes from rs data                
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
                    Si.mvc          = Bj.mvc;
                    Si.FM           = Bj.FM;
                    Si.ARAT         = Bj.ARAT;
                    Si.lesion_PrCG  = SSi.lesion_PrCG(1);
                    Si.lesion_PrCG_W= SSi.lesion_PrCG_W(1);
                else
                    Si.ensOverall   = nan;
                    Si.mmOverall    = nan;
                    Si.mvcNorm      = nan;
                    Si.mvc          = nan;
                    Si.FM           = nan;
                    Si.ARAT         = nan;
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

                % clean out regions with no signal
                tsClean = ts;
                flag    = 0;
                for tsi = 1:length(ts)
                    tsClean_i       = tsClean{tsi};
                    idx             = std(tsClean_i,[],1);
                    idx             = idx~=0 & ~isnan(idx);
                    tsClean{tsi}    = tsClean_i(:,idx);
                    if isempty(tsClean{tsi})
                        flag = 1;
                    end;
                end
                if flag==1
                    continue;
                end;
                

                % depending on whether it is a control or a patient, set
                % the relevant reference region
                %   - left M1 for controls
                %   - lesioned M1 for patients
                refTS   = []; % reference time-series (lesioned M1 for patients, left hem for controls)
                compTS  = [];
                if (Si.control==1)
                    refTS   = tsClean(1:2);
                    compTS  = tsClean(3:4);
                else
                    if (Si.lesionSide==1)
                        refTS   = tsClean(1:2);
                        compTS  = tsClean(3:4);
                    else
                        refTS   = tsClean(3:4);
                        compTS  = tsClean(1:2);
                    end;
                end;
                
                % avg correlations of reference with itself
                cRef = fisherz(corr(refTS{2}));
                cRef(isinf(cRef)) = nan;
                cRef = fisherinv(nanmean(cRef,2));
                
                % avg correlations of reference with whole brain 
                wb  = [compTS{2} compTS{1} refTS{1}];
                cComp = fisherz(corr(refTS{2},wb));
                cComp(isinf(cComp)) = nan;
                cComp = fisherinv(nanmean(cComp,1))';
                
                % estimate 95% percentile of correlations
                pr = prctile([cRef;cComp],95);
                
                % prc voxels in reference M1
                voxRef = sum(cRef>pr)/length(cRef);
                
                % prc voxels in contra M1
                vox     = cComp(1:size(compTS{2},2));
                voxComp = sum(vox>pr)/length(vox);
                
                % laterality index
                LI = voxRef - voxComp;
                                
                Si.voxRef   = voxRef;
                Si.voxComp  = voxComp;
                Si.LI       = LI;
                
                S           = addstruct(S,Si);
            end;
        end;
        varargout = {S};
        save(fullfile(ppDir,'park2011_preprocess.mat'),'-struct','S');   
        
    case 'REP_LI_Park2011'
        pre = {'preproc_Bold_Rest'};
        lesionType = {'none','cortical','subcortical','mixed'}; %1,2,3,4
        roi  = combnk(1:16,2); % all pairs
        inclROI = [2 10];
        
        % 0. check which subjects have the raw and preprocessed files
        D = rs_preprocess_ana('check_data_exists');
        D = getrow(D,D.existPPImg & D.existRawImg);

        % remove bad control session
        idx = ismember(D.subj_name,'CUP_1002') & D.week==24;
        D   = getrow(D,~idx);        

%         D   = dload(fullfile(ppDir,'patient_list.txt'));
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
                R.R = R.R(inclROI);
                fDir= fullfile(rsDir,D.subj_name{i},D.Week{i});
                v   = spm_vol(fullfile(fDir,fname));       % load volumes from rs data
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
                    Si.mvc          = Bj.mvc;
                    Si.FM           = Bj.FM;
                    Si.ARAT         = Bj.ARAT;
                    Si.lesion_PrCG  = SSi.lesion_PrCG(1);
                    Si.lesion_PrCG_W= SSi.lesion_PrCG_W(1);
                else
                    Si.ensOverall   = nan;
                    Si.mmOverall    = nan;
                    Si.mvcNorm      = nan;
                    Si.mvc          = nan;
                    Si.FM           = nan;
                    Si.ARAT         = nan;
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

                % estimate left and right combined M1/S1 regions
                ts1 = [ts{1},ts{2}];    % left hemisphere (M1/S1)
                ts2 = [ts{3},ts{4}];    % right hemisphere (M1/S1)
                
                % removing all voxels with no signal
                ts1 = ts1(:,all(ts1,1));    
                ts2 = ts2(:,all(ts2,1));    

                refTS   = []; % reference time-series (lesioned hem for patients, left hem for controls)
                compTS  = [];
                if (Si.control==1)
                    refTS   = ts1;
                    compTS  = ts2;
                else
                    if (Si.lesionSide==1)
                        refTS   = ts1;
                        compTS  = ts2;
                    else
                        refTS   = ts2;
                        compTS  = ts1;
                    end;
                end;
                
                % estimate pair-wise correlations
                %   - fisherz transform (distribution exists around -3 and 3)
                %   - remove inf
                r       = corr(refTS);
                idx     = any(~isnan(r),1);
                refFZ   = fisherz(1-squareform(1-r(idx,idx)));
                refFZ   = refFZ(refFZ>-3 & refFZ<3);
                
                r       = corr(refTS,compTS);
                compFZ  = fisherz(r(~isnan(r)));
                compFZ  = compFZ(:);
                compFZ  = compFZ(compFZ>-3 & compFZ<3);
                
                Si.compFZ   = mean(compFZ);
                Si.refFZ    = mean(refFZ);
                Si.refCon   = mean(Si.compFZ)./mean(Si.refFZ);
                
                S           = addstruct(S,Si);
            end;
        end;
        varargout = {S};
        save(fullfile(ppDir,'park2013_preprocess.mat'),'-struct','S'); 
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
    case 'REL_withinSubj_week'
       
        fName = {'nc_interhem','nc_intrahem_lesioned','nc_intrahem_nonlesioned','nc_getThemall'};

for i=1:length(fName)
        D = load(fullfile([ppDir],sprintf('%s.mat',fName{i})));

        weeks = {[1 4],[1 12],[1 24],[1 52]};
        %weeks = {[1 4],[4 12],[12 24],[24 52]};

    S = [];
    for sn=unique(D.subj_name)'
    Ds = getrow(D,ismember(D.subj_name,sn));
    
        for w=1:length(weeks)
        weekOfInterest = weeks{w};
        
        Dw1 = getrow(Ds,Ds.week==weekOfInterest(1));
        Dw2 = getrow(Ds,Ds.week==weekOfInterest(2));
        
            if (isempty(Dw1.subj_name) || isempty(Dw2.subj_name))
%             continue;
%         end;
            Si                  = getrow(Ds,1);
            Si.r                = NaN;
            Si.fz_r             = NaN;
            Si.weekOfInterest   = weekOfInterest;
            Si.weekLabel        = w;
            S                   = addstruct(S,Si); 
            else
            Si                  = getrow(Ds,1);
            Si.r                = corr(Dw1.M',Dw2.M');
            Si.fz_r             = fisherz(Si.r);
            Si.weekOfInterest   = weekOfInterest;
            Si.weekLabel        = w;
            S                   = addstruct(S,Si);
            end;
        end;
    end;
        
        varargout = {S};
        
   
%          subplot(2,2,i)
%          s   = style_sheet(sty_2grp_cp,'leg',{'patient','control'},'leglocation','northoutside','errorcolor',sty_2grp_cp);
%          
%           barplot(S.weekLabel,S.r,'plotfcn','fisherinv(nanmean(x))','split',S.control,'CAT',s(1).CAT,s(1).PLOT{:});
%           set_graphics(gcf,'ylabel',{'correlation'}, 'xlabel',{'weeks'},'fontsize_all',14);
%           title(fName(i));          
%           % [~,~,sn] = unique(S.subj_name);
%           % anovaMixed(S.r,sn,'between',S.control,{'grp'});
          
        
          subplot(2,2,i)
          s   = style_sheet(sty_2grp_cp,'leg',{'patient','control'},'leglocation','northoutside','errorcolor',sty_2grp_cp);
         
          barplot(S.weekLabel,S.r,'plotfcn','fisherinv(nanmean(x))','split',S.control,'CAT',s(1).CAT,s(1).PLOT{:});
          set_graphics(gcf,'ylabel',{'correlation'}, 'xlabel',{'weeks'},'fontsize_all',14);
          title(fName(i));
        
          % test differences between weeLabel
        SS = getrow(S,S.control==0);  
        x = pivottable(SS.subj_name,SS.weekLabel,SS.fz_r,'nanmean');
        ttest(x(:,1),x(:,2),2,'paired');
        ttest(x(:,1),x(:,3),2,'paired');
        
        P = getrow(S,S.control==0);
             m=nanmean(P.r); 
             SE=nanstd(P.r)./sqrt(sum(~isnan(P.r),1));
             fprintf('Inter-subject patient r = %2.3f (%2.3f - %2.3f)\n',(m),(m-1.96*SE),(m+1.96*SE));
        C = getrow(S,S.control==1);
             m=nanmean(C.r); 
             SE=nanstd(C.r)./sqrt(sum(~isnan(C.r),1));
             fprintf('Inter-subject control r = %2.3f (%2.3f - %2.3f)\n',(m),(m-1.96*SE),(m+1.96*SE));
             
         name = char(fName(i));
         SS =     tapply(S,{'subj_name','weekLabel','control'},{'fz_r','nanmean'},{'r','nanmean'});
         save(fullfile(ppDir,[name,'_all.mat']),'-struct','SS');   
         dsave(fullfile(statsDir,[name,'_all.dat']),SS);
                   
end;
   keyboard; 
   
   
   case 'REL_withinSubj_week_EUC_control'
       
        fName = {'nc_interhem','nc_intrahem_lesioned','nc_intrahem_nonlesioned','nc_getThemall'};

    for i=1:length(fName)
        D = load(fullfile([ppDir],sprintf('%s.mat',fName{i})));

        weeks = {[1 4],[1 12],[1 24],[1 52]};
        %weeks = {[1 4],[4 12],[12 24],[24 52]};

    S = [];
    for sn=unique(D.subj_name)'
    Ds = getrow(D,ismember(D.subj_name,sn));
    
        for w=1:length(weeks)
        weekOfInterest = weeks{w};
        
        Dw1 = getrow(Ds,Ds.week==weekOfInterest(1));
        Dw2 = getrow(Ds,Ds.week==weekOfInterest(2));
        
            if (isempty(Dw1.subj_name) || isempty(Dw2.subj_name))
%             continue;
%         end;
            Si                  = getrow(Ds,1);
            Si.r                = NaN;
            Si.fz_r             = NaN;
            Si.weekOfInterest   = weekOfInterest;
            Si.weekLabel        = w;
            S                   = addstruct(S,Si); 
            else
            Si                  = getrow(Ds,1);
            Si.r                = corr(Dw1.M',Dw2.M');
            Si.fz_r             = fisherz(Si.r);
            Si.weekOfInterest   = weekOfInterest;
            Si.weekLabel        = w;
            S                   = addstruct(S,Si);
            end;
        end;
    end;
        
        varargout = {S};
        
   
%          subplot(2,2,i)
%          s   = style_sheet(sty_2grp_cp,'leg',{'patient','control'},'leglocation','northoutside','errorcolor',sty_2grp_cp);
%          
%           barplot(S.weekLabel,S.r,'plotfcn','fisherinv(nanmean(x))','split',S.control,'CAT',s(1).CAT,s(1).PLOT{:});
%           set_graphics(gcf,'ylabel',{'correlation'}, 'xlabel',{'weeks'},'fontsize_all',14);
%           title(fName(i));          
%           % [~,~,sn] = unique(S.subj_name);
%           % anovaMixed(S.r,sn,'between',S.control,{'grp'});
          
        
          subplot(2,2,i)
          s   = style_sheet(sty_2grp_cp,'leg',{'patient','control'},'leglocation','northoutside','errorcolor',sty_2grp_cp);
         
          barplot(S.weekLabel,S.r,'plotfcn','fisherinv(nanmean(x))','split',S.control,'CAT',s(1).CAT,s(1).PLOT{:});
          set_graphics(gcf,'ylabel',{'correlation'}, 'xlabel',{'weeks'},'fontsize_all',14);
          title(fName(i));
        
          % test differences between weeLabel
        SS = getrow(S,S.control==0);  
        x = pivottable(SS.subj_name,SS.weekLabel,SS.fz_r,'nanmean');
        ttest(x(:,1),x(:,2),2,'paired');
        ttest(x(:,1),x(:,3),2,'paired');
        
        P = getrow(S,S.control==0);
             m=nanmean(P.r); 
             SE=nanstd(P.r)./sqrt(sum(~isnan(P.r),1));
             fprintf('Inter-subject patient r = %2.3f (%2.3f - %2.3f)\n',(m),(m-1.96*SE),(m+1.96*SE));
        C = getrow(S,S.control==1);
             m=nanmean(C.r); 
             SE=nanstd(C.r)./sqrt(sum(~isnan(C.r),1));
             fprintf('Inter-subject control r = %2.3f (%2.3f - %2.3f)\n',(m),(m-1.96*SE),(m+1.96*SE));
             
         name = char(fName(i));
         SS =     tapply(S,{'subj_name','weekLabel','control'},{'fz_r','nanmean'},{'r','nanmean'});
         save(fullfile(ppDir,[name,'_all.mat']),'-struct','SS');   
         dsave(fullfile(statsDir,[name,'_all.dat']),SS);
                   
end;
   keyboard;     
    case 'REL_withinAll'                % calculate within subject correlations (early/late after stroke)
        fName = {'nc_interhem','nc_intrahem_lesioned','nc_intrahem_nonlesioned','nc_getThemall'};
        
        S = [];
        for i=1:length(fName)
            D = load(fullfile([ppDir],sprintf('%s.mat',fName{i})));
            D = rs_imana_110118('REL_withinSubj',D);
            
            Si      = D;
            Si.type = zeros(length(D.subj_name),1)+i;
            S       = addstruct(S,Si);
            
            disp(fName{i});
            x = pivottable(D.subj_name,D.control,D.fz_r,'nanmean');
            ttest(x(:,1),x(:,2),2,'independent');
        end;
        
        
        h = make_figure;
        s   = style_sheet(sty_2grp_cp,'leg',{'patient','control'},'leglocation','northoutside','errorcolor',sty_2grp_cp);
        
        barplot(S.type,S.fz_r,'plotfcn','fisherinv(mean(x))','split',S.control,'CAT',s(1).CAT,s(1).PLOT{:});
        set_graphics(gcf,'ylabel',{'correlation acute/chronic'}, 'xlabel',{'type'},'fontsize_all',14);
              
        [~,~,sn] = unique(S.subj_name);
        anovaMixed(S.fz_r,sn,'between',S.control,{'grp'},'within',S.type,{'type'});
        keyboard;
    case 'REL_fingerprint'              % check to see if a finger print exists for the participants
        fName = {'nc_interhem','nc_intrahem_lesioned','nc_intrahem_nonlesioned'};
        
        S = [];
        for i=1:length(fName)
            D = load(fullfile(rootDir,'/new',sprintf('%s.mat',fName{i})));
            D = getrow(D,D.lesionType ~=2);
            D = getrow(D,D.lesionType ~=4);
            D = rs_imana_110118('REL_withinSubj',D);
            
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
        
        h   = make_figure;
        s   = style_sheet(sty_multiple,'leg',{'within','across early','across late'},'leglocation','northoutside','errorcolor',sty_multiple);
        
        subplot(121);
        barplot(T.type,T.fz,'split',T.fp,'subset',T.control==1,'CAT',s(1).CAT,s(1).PLOT{:});
        title('control');
        subplot(122);
        barplot(T.type,T.fz,'split',T.fp,'subset',T.control==0,'CAT',s(1).CAT,s(1).PLOT{:});
        title('patient');
        set_graphics(gcf,'ylabel',{'correlation'}, 'xlabel',{'type'},'fontsize_all',14);
        
        
        [~,~,T.sn] = unique(T.subj_name);
        T1 = getrow(T,T.control==1);
        disp('controls');
        anovaMixed(T1.fz,T1.sn,'within',[T1.type T1.fp],{'fingerprint','type'});
        
        T2 = getrow(T,T.control==0);
        disp('patients');
        anovaMixed(T2.fz,T2.sn,'within',[T2.type T2.fp],{'fingerprint','type'});
   keyboard;
    case 'PP_patientSplit'              % get subset split of patients with clinical characteristics
        Dp = varargin{1};
        vararginoptions({varargin{2:end}},{'lesionType'});
        
        % 0. Avg controls (all) and patients (per week and according to severity)
        C           = tapply(D,{'subj_name'},{'mat','mean(x,1)','subset',D.control==1});

        D.mvcNorm   = -D.mvcNorm;                       % for some reason mvcNorm is negative (to be fixed)
        [x,sn]      = pivottable(D.subj_name,D.week,D.mvcNorm,'mean(x)','subset',D.control==0);
        x           = nanmean(x(:,1:2),2);

    case 'FIG_interhemSomat'                            % somatotopic ordering of correlations between hemispheres    
        S = load(fullfile(ppDir,'nc_interhem.mat'));
        S = getrow(S,S.lesionType ~=2);
        S = getrow(S,S.lesionType ~=4);
        
        hom = 1:6:25;                       % homologous interhemispheric areas
        het = find(~ismember(1:25,hom));    % heterologous interhemispheric areas
        
        % get required roi pairs
        S.hom = nanmean(S.M(:,hom),2);
        S.het = nanmean(S.M(:,het),2);
        
        % make plot for hom/het pairs
        s   = style_sheet(sty_2grp_ac,'leg',{'hom','het'},'leglocation','northoutside','errorcolor',sty_2grp_ac);
    
        h = make_figure;
        subplot(121);
        lineplot(S.week,[S.hom S.het],'plotfcn','fisherinv(nanmean(x))','subset',S.control==1,'CAT',s(1).CAT,s(1).PLOT{:});
        title('controls');
        ylim([0.2 1.1]);
        subplot(122);
        lineplot(S.week,[S.hom S.het],'plotfcn','fisherinv(nanmean(x))','subset',S.control==0,'CAT',s(1).CAT,s(1).PLOT{:});
        title('patients');
        ylim([0.2 1.1]);
        set_graphics(gcf,'ylabel',{'FC (fisherz)'}, 'xlabel',{'weeks'},'fontsize_all',14);
        
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
        dsave(fullfile(statsDir,[what,'.dat']),D);
    case 'FIG_inter_splithalf_reliability'              % split half reliability of patterns within each individual/session
         S = load(fullfile(ppDir,'nc_interhem_splithalf.mat'));        
         S = getrow(S,S.lesionType ~=2);
         S = getrow(S,S.lesionType ~=4);
%         
        s   = style_sheet(sty_2grp_cp,'errorcolor',sty_2grp_cp);
        barplot([S.week],S.sh_r,'split',S.control,'plotfcn','fisherinv(nanmean(x))','CAT',s(1).CAT,s(1).PLOT{:});
        set_graphics(gcf,'ylabel',{'split-half correlation'},'fontsize_all',14);
        title('interhemispheric')
        keyboard;
        
        C = getrow(S,S.control==1);
        P = getrow(S,S.control==0);
            m=mean(C.sh_r); 
            SE=stderr(C.sh_r); 
            fprintf('Inter-subject r = %2.3f (%2.3f - %2.3f)\n',(m),(m-1.96*SE),(m+1.96*SE));
        
            m=mean(P.sh_r); 
            SE=stderr(P.sh_r); 
            fprintf('Inter-subject r = %2.3f (%2.3f - %2.3f)\n',(m),(m-1.96*SE),(m+1.96*SE));

        % save for statistics
        D = tapply(S,{'subj_name','week','control'},{'sh_fz_r','nanmean'});
        dsave(fullfile(statsDir,[what,'.dat']),D);
    case 'FIG_intraLesioned_splithalf_reliability'      % split half reliability of patterns within each individual/session
        S = load(fullfile(ppDir,'nc_intraHem_lesioned_splithalf.mat'));        
        S = getrow(S,S.lesionType ~=2);
        S = getrow(S,S.lesionType ~=4);
        
        s   = style_sheet(sty_2grp_cp,'errorcolor',sty_2grp_cp);
        barplot([S.week],S.sh_r,'split',S.control,'plotfcn','fisherinv(nanmean(x))','CAT',s(1).CAT,s(1).PLOT{:});
        set_graphics(gcf,'ylabel',{'split-half correlations)'},'fontsize_all',14);
        keyboard
        pivottable([],S.control,S.sh_fz_r,'nanmean(fisherinv(nanmean(x)))');
        
        C = getrow(S,S.control==1);
        P = getrow(S,S.control==0);
            m=mean(C.sh_r); 
            SE=stderr(C.sh_r); 
            fprintf('Inter-subject r = %2.3f (%2.3f - %2.3f)\n',(m),(m-1.96*SE),(m+1.96*SE));
        
            m=mean(P.sh_r); 
            SE=stderr(P.sh_r); 
            fprintf('Inter-subject r = %2.3f (%2.3f - %2.3f)\n',(m),(m-1.96*SE),(m+1.96*SE));

        % save for statistics
        D = tapply(S,{'subj_name','week','control'},{'sh_fz_r','nanmean'});
        dsave(fullfile(statsDir,[what,'.dat']),D);        
    case 'FIG_intraNonlesioned_splithalf_reliability'   % split half reliability of patterns within each individual/session
        S = load(fullfile(ppDir,'nc_intraHem_nonlesioned_splithalf.mat'));        
        S = getrow(S,S.lesionType ~=2);
        S = getrow(S,S.lesionType ~=4);
        
        s   = style_sheet(sty_2grp_cp,'errorcolor',sty_2grp_cp);
        set_graphics(gcf,'ylabel',{'split-half correlation'},'fontsize_all',14);
        barplot([S.week],S.sh_r,'split',S.control,'plotfcn','fisherinv(nanmean(x))','CAT',s(1).CAT,s(1).PLOT{:});
        set_graphics(gcf,'ylabel',{'FC (fisherz)'}, 'xlabel',{'weeks'},'fontsize_all',14);
        keyboard
        pivottable([],S.control,S.sh_fz_r,'nanmean(fisherinv(nanmean(x)))');
        
        C = getrow(S,S.control==1);
        P = getrow(S,S.control==0);
            m=mean(C.sh_r); 
            SE=stderr(C.sh_r); 
            fprintf('Inter-subject r = %2.3f (%2.3f - %2.3f)\n',(m),(m-1.96*SE),(m+1.96*SE));
        
            m=mean(P.sh_r); 
            SE=stderr(P.sh_r); 
            fprintf('Inter-subject r = %2.3f (%2.3f - %2.3f)\n',(m),(m-1.96*SE),(m+1.96*SE));

        % save for statistics
        D = tapply(S,{'subj_name','week','control'},{'sh_fz_r','nanmean'});
        dsave(fullfile(statsDir,[what,'.dat']),D);                
    case 'FIG_splithalf_all'
         S = load(fullfile(ppDir,'nc_getThemall_splithalf.mat'));        
         S = getrow(S,S.lesionType ~=2);
         S = getrow(S,S.lesionType ~=4);
%         
        s   = style_sheet(sty_2grp_cp,'errorcolor',sty_2grp_cp);
        barplot([S.week],S.sh_r,'split',S.control,'plotfcn','fisherinv(nanmean(x))','CAT',s(1).CAT,s(1).PLOT{:});
        set_graphics(gcf,'ylabel',{'split-half correlation'},'fontsize_all',14);
        
        C = getrow(S,S.control==1);
        P = getrow(S,S.control==0);
            m=mean(C.sh_r); 
            SE=stderr(C.sh_r); 
            fprintf('Inter-subject r = %2.3f (%2.3f - %2.3f)\n',(m),(m-1.96*SE),(m+1.96*SE));
        
            m=mean(P.sh_r); 
            SE=stderr(P.sh_r); 
            fprintf('Inter-subject r = %2.3f (%2.3f - %2.3f)\n',(m),(m-1.96*SE),(m+1.96*SE));

        % save for statistics
        D = tapply(S,{'subj_name','week','control'},{'sh_fz_r','nanmean'});
        dsave(fullfile(statsDir,[what,'.dat']),D);
    case 'FIG_all'    
 
        subplot(221);
        rs_imana_110118('FIG_inter_splithalf_reliability');
        ylim([0 1]);
        ylabel('Fisher Z')
        xlabel('weeks')
        title('Interhemispheric')
        subplot(222);
        rs_imana_110118('FIG_intraLesioned_splithalf_reliability');
        ylim([0 1]);
        ylabel('Fisher Z')
        xlabel('weeks')
        title('Intrahemispheric Les')
        subplot(223);
        rs_imana_110118('FIG_intraNonlesioned_splithalf_reliability');
        ylim([0 1]); 
        ylabel('Fisher Z')
        xlabel('weeks')
        title('Intrahemispheric NonLes')
        subplot(224);
        rs_imana_110118('FIG_splithalf_all');
        ylim([0 1]); 
        ylabel('Fisher Z')
        xlabel('weeks')
        title('All connections')
    case 'FIG_RelCon'
     D   = load(fullfile(rootDir,'golestani2013_preprocess.mat'));   
 
      % 1. get controls 
        Dc      = getrow(D,D.control==1);
        [x,sn]  = pivottable(Dc.subj_name,[],Dc.week,'length(unique(x))');
        Dc      = getrow(Dc,ismember(Dc.subj_name,sn(x>1)));
        
        Dp      = getrow(D,D.lesionType==3);
        [x,sn]  = pivottable(Dp.subj_name,[],Dp.week,'length(unique(x))');
        Dp      = getrow(Dp,ismember(Dp.subj_name,sn(x>1)));
        
        T = addstruct(Dc,Dp);
        varargout = {T};
     
     h = make_figure;
     s   = style_sheet(sty_2grp_ac,'leg',{'','het'},'leglocation','northoutside','errorcolor',sty_2grp_ac);
    
     title('RelCon for interhemispheric Sensorimotor cortex (S1M1-S1M1)')
     lineplot(T.week,T.refCon,'plotfcn','nanmean','errorfcn','stderr','split',T.control,'CAT',s(1).CAT,s(1).PLOT{:});
     ylabel('RelCon')
     xlabel('weeks')
     
     % correlation between RelCon and FM
     % needs to be done for hand strength and ARAT too
     
     Tp = getrow(T,T.control==0);
     
     T1     = getrow(Tp,Tp.week==1);
     T4     = getrow(Tp,Tp.week==4);
     T12    = getrow(Tp,Tp.week==12);
     T24    = getrow(Tp,Tp.week==24);
     T52    = getrow(Tp,Tp.week==52);
     
     [r1,p1] = corr(T1.refFZ,T1.FM,'row','complete');
     [r2,p2] = corr(T4.refFZ,T4.FM,'row','complete');
     [r3,p3] = corr(T12.refFZ,T12.FM,'row','complete');
     [r4,p4] = corr(T24.refFZ,T24.FM,'row','complete');
     [r5,p5] = corr(T52.refFZ,T52.FM,'row','complete');
     
        DD =     tapply(T,{'subj_name','week','control'},{'refCon'});
        dsave(fullfile(statsDir,'RelCon_stat.dat'),DD); 
        
        keyboard;
        
    case 'FIG_Park'
    D   = load(fullfile(ppDir,'park2011_preprocess.mat'));  
              
         Dc      = getrow(D,D.control==1);
        [x,sn]  = pivottable(Dc.subj_name,[],Dc.week,'length(unique(x))');
        Dc      = getrow(Dc,ismember(Dc.subj_name,sn(x>1)));
        
        Dp      = getrow(D,D.lesionType==3);
        [x,sn]  = pivottable(Dp.subj_name,[],Dp.week,'length(unique(x))');
        Dp      = getrow(Dp,ismember(Dp.subj_name,sn(x>1)));
        
        T = addstruct(Dc,Dp);
        varargout = {T};
        
     h = make_figure;
     s   = style_sheet(sty_2grp_ac,'leg',{'patients','controls'},'leglocation','northoutside','errorcolor',sty_2grp_ac);
    
     title('Laterality_Index (M1les-M1nonles)')
     lineplot(T.week,T.LI,'plotfcn','nanmean','errorfcn','stderr','split',T.control,'CAT',s(1).CAT,s(1).PLOT{:});
     ylabel('LI')
     xlabel('weeks')
     
     DD =     tapply(T,{'subj_name','week','control'},{'LI'});
     dsave(fullfile(statsDir,'LI_stat.dat'),DD); 
        
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
        D       = rs_imana_110118('PP_studySubset'); 
        D       = rs_imana_110118('PP_removePair',D);
    
    % only take subcorticals for now
      D   = getrow(D,(D.lesionType~=2)); %no cortical
      D   = getrow(D,(D.lesionType~=4)); %no mixed
    
%     or define threshold
%     D   = getrow(D,D.lesion_PrCG<=0.05);
    
    
%  Roi pairs
     
    M1_M1    = getrow(D,D.rType==23);        %[2 10] 
%      Pmd_Pmd  = getrow(D,D.rType==37);        %[3 11]
%      Pmv_Pmv  = getrow(D,D.rType==50);        %[4 12]
%      SMA_SMA  = getrow(D,D.rType==62);        %[5 13]
%      V1_V1    = getrow(D,D.rType==73);        %[6 14]

    M1 = tapply(M1_M1,{'subj_name','week','control'},{'fz_r','nanmean(x,1)'});
    dsave(fullfile(statsDir,'M1_timecourse.dat'),M1);
    
    
    s   = style_sheet(sty_2grp_ac,'leg',{'patients','control'},'leglocation','northoutside');
     
    %interhemispheric connectivity 
    h = make_figure;
%     subplot(3,2,1)
    title('M1-M1 subcortical')
    %lineplot(M1_M1.week,M1_M1.fz_r,'plotfcn','nanmean','errorfcn','stderr','split',M1_M1.control,'CAT',s(1).CAT,s(1).PLOT{:});
    lineplot(M1_M1.week,M1_M1.r,'plotfcn','nanmean','errorfcn','stderr','split',M1_M1.control,'CAT',s(1).CAT,s(1).PLOT{:});
    ylabel('correlation index')
    xlabel('weeks')
%     subplot(3,2,2)
%     title('Pmd-Pmd subcortical')
%     %lineplot(Pmd_Pmd.week,Pmd_Pmd.fz_r,'plotfcn','nanmean','errorfcn','stderr','split',Pmd_Pmd.control,'CAT',s(1).CAT,s(1).PLOT{:});
%     lineplot(Pmd_Pmd.week,Pmd_Pmd.r,'plotfcn','nanmean','errorfcn','stderr','split',Pmd_Pmd.control,'CAT',s(1).CAT,s(1).PLOT{:});
%     ylabel('correlation index')
%     xlabel('weeks')
%     subplot(3,2,3)
%     title('Pmv-Pmv subcortical')
%     %lineplot(Pmv_Pmv.week,Pmv_Pmv.fz_r,'plotfcn','nanmean','errorfcn','stderr','split',Pmv_Pmv.control,'CAT',s(1).CAT,s(1).PLOT{:});
%     lineplot(Pmv_Pmv.week,Pmv_Pmv.r,'plotfcn','nanmean','errorfcn','stderr','split',Pmv_Pmv.control,'CAT',s(1).CAT,s(1).PLOT{:});
%     ylabel('correlation index')
%     xlabel('weeks')
%     subplot(3,2,4)
%     title('SMA-SMA subcortical')
%     %lineplot(SMA_SMA.week,SMA_SMA.fz_r,'plotfcn','nanmean','errorfcn','stderr','split',SMA_SMA.control,'CAT',s(1).CAT,s(1).PLOT{:});
%     lineplot(SMA_SMA.week,SMA_SMA.r,'plotfcn','nanmean','errorfcn','stderr','split',SMA_SMA.control,'CAT',s(1).CAT,s(1).PLOT{:});
%     ylabel('correlation index')
%     xlabel('weeks') 
%     subplot(3,2,5)
%     title('V1-V1 subcortical')
%     %lineplot(V1_V1.week,V1_V1.fz_r,'plotfcn','nanmean','errorfcn','stderr','split',V1_V1.control,'CAT',s(1).CAT,s(1).PLOT{:});
%     lineplot(V1_V1.week,V1_V1.r,'plotfcn','nanmean','errorfcn','stderr','split',V1_V1.control,'CAT',s(1).CAT,s(1).PLOT{:});
%     ylabel('correlation index')
%     xlabel('weeks') 
    
    %% interhemispheric individual patients
    p   = style_sheet(sty_multiple,'leg','auto','leglocation','northoutside');
    
    h = make_figure;
%     subplot(3,2,1)
    title('M1-M1 subcortical')
    %lineplot(M1_M1.week,M1_M1.fz_r,'plotfcn','nanmean','errorfcn','stderr','split',M1_M1.subj_name,'subset',M1_M1.control==0);
    lineplot(M1_M1.week,M1_M1.r,'plotfcn','nanmean','errorfcn','stderr','split',M1_M1.subj_name,'subset',M1_M1.control==0);
    ylabel('correlation index')
    xlabel('weeks')
%     subplot(3,2,2)
%     title('Pmd-Pmd subcortical')
%     lineplot(Pmd_Pmd.week,Pmd_Pmd.fz_r,'plotfcn','nanmean','errorfcn','stderr','split',Pmd_Pmd.subj_name,'subset',Pmd_Pmd.control==0);
%     %lineplot(Pmd_Pmd.week,Pmd_Pmd.r,'plotfcn','nanmean','errorfcn','stderr','split',Pmd_Pmd.subj_name,'subset',Pmd_Pmd.control==0);
%     ylabel('correlation index')
%     xlabel('weeks')
%     subplot(3,2,3)
%     title('Pmv-Pmv subcortical')
%     %lineplot(Pmv_Pmv.week,Pmv_Pmv.fz_r,'plotfcn','nanmean','errorfcn','stderr','split',Pmv_Pmv.subj_name,'subset',Pmv_Pmv.control==0);
%     lineplot(Pmv_Pmv.week,Pmv_Pmv.r,'plotfcn','nanmean','errorfcn','stderr','split',Pmv_Pmv.subj_name,'subset',Pmv_Pmv.control==0);
%     ylabel('correlation index')
%     xlabel('weeks')
%     subplot(3,2,4)
%     title('SMA-SMA subcortical')
%     %lineplot(SMA_SMA.week,SMA_SMA.fz_r,'plotfcn','nanmean','errorfcn','stderr','split',SMA_SMA.subj_name,'subset',SMA_SMA.control==0);
%     lineplot(SMA_SMA.week,SMA_SMA.r,'plotfcn','nanmean','errorfcn','stderr','split',SMA_SMA.subj_name,'subset',SMA_SMA.control==0);
%     ylabel('correlation index')
%     xlabel('weeks')
%     subplot(3,2,5)
%     title('V1-V1 subcortical')
%     %lineplot(V1_V1.week,V1_V1.fz_r,'plotfcn','nanmean','errorfcn','stderr','split',V1_V1.subj_name,'subset',V1_V1.control==0);
%     lineplot(V1_V1.week,V1_V1.r,'plotfcn','nanmean','errorfcn','stderr','split',V1_V1.subj_name,'subset',V1_V1.control==0);
%     ylabel('correlation index')
%     xlabel('weeks')
    
% %% intracortical connectivity/ lesioned/dominant
%     
%      h = make_figure;
% %      subplot(1,4,1)
%      title('M1-SMA subcortical')
%      M1_SMA = getrow(D, D.regType(:,1)==2 & D.regType(:,2)==5);
%      lineplot(M1_SMA.week,M1_SMA.r,'plotfcn','nanmean','errorfcn','stderr','split',M1_SMA.control,'CAT',s(1).CAT,s(1).PLOT{:}); 
%      %individuals
%      %lineplot(M1_SMA.week,M1_SMA.r,'plotfcn','nanmean','errorfcn','stderr','split',M1_SMA.subj_name,'subset',M1_SMA.control==0); 
%      ylabel('correlation index')
%      xlabel('weeks')
%      subplot(1,4,2)
%      title('M1-Pmd subcortical')
%      M1_Pmd = getrow(D, D.regType(:,1)==2 & D.regType(:,2)==3);
%      lineplot(M1_Pmd.week,M1_Pmd.r,'plotfcn','nanmean','errorfcn','stderr','split',M1_Pmd.control,'CAT',s(1).CAT,s(1).PLOT{:});
%      %individuals
%      %lineplot(M1_Pmd.week,M1_Pmd.r,'plotfcn','nanmean','errorfcn','stderr','split',M1_Pmd.subj_name,'subset',M1_Pmd.control==0);
%      ylabel('correlation index')
%      xlabel('weeks')
%      subplot(1,4,3)
%      title('M1-Pmv subcortical')
%      M1_Pmv = getrow(D, D.regType(:,1)==2 & D.regType(:,2)==4);
%      lineplot(M1_Pmv.week,M1_Pmv.r,'plotfcn','nanmean','errorfcn','stderr','split',M1_Pmv.control,'CAT',s(1).CAT,s(1).PLOT{:});
%      %individuals
%      %lineplot(M1_Pmv.week,M1_Pmv.r,'plotfcn','nanmean','errorfcn','stderr','split',M1_Pmv.subj_name,'subset',M1_Pmd.control==0); 
%      ylabel('correlation index')
%      xlabel('weeks')
%      subplot(1,4,4)
%      title('M1-S1 subcortical')
%      M1_S1 = getrow(D, D.regType(:,1)==1 & D.regType(:,2)==2);
%      lineplot(M1_S1.week,M1_S1.r,'plotfcn','nanmean','errorfcn','stderr','split',M1_S1.control,'CAT',s(1).CAT,s(1).PLOT{:});
%      %individuals
%      %lineplot(M1_S1.week,M1_S1.r,'plotfcn','nanmean','errorfcn','stderr','split',M1_S1.subj_name,'subset',M1_S1.control==0);
%      ylabel('correlation index')
%      xlabel('weeks')
%      
% %% intracortical connectivity/ NON lesioned/dominant
%      T = getrow(D, D.lesionSide==D.hl);     
%      h = make_figure;
%      subplot(1,4,1)
%      title('M1-SMA subcortical')
%      M1_SMA = getrow(T, T.regType(:,1)==2 & T.regType(:,2)==5);
%      lineplot(M1_SMA.week,M1_SMA.r,'plotfcn','nanmean','errorfcn','stderr','split',M1_SMA.control); ylabel('correlation index')
%      xlabel('weeks')
%      subplot(1,4,2)
%      title('M1-Pmd subcortical')
%      M1_Pmd = getrow(T, T.regType(:,1)==2 & T.regType(:,2)==3);
%      lineplot(M1_Pmd.week,M1_Pmd.r,'plotfcn','nanmean','errorfcn','stderr','split',M1_Pmd.control,'CAT',s(1).CAT,s(1).PLOT{:});
%      ylabel('correlation index')
%      xlabel('weeks')
%      subplot(1,4,3)
%      title('M1-Pmv subcortical')
%      M1_Pmv = getrow(T, T.regType(:,1)==2 & T.regType(:,2)==4);
%      lineplot(Pmv.week,Pmv.r,'plotfcn','nanmean','errorfcn','stderr','split',M1.control,'CAT',s(1).CAT,s(1).PLOT{:});
%      ylabel('correlation index')
%      xlabel('weeks')
%      subplot(1,4,4)
%      title('M1-S1 subcortical')
%      M1_S1 = getrow(T, T.regType(:,1)==2 & T.regType(:,2)==1);
%      lineplot(M1_S1.week,M1_S1.r,'plotfcn','nanmean','errorfcn','stderr','split',M1_S1.control,'CAT',s(1).CAT,s(1).PLOT{:});
%      ylabel('correlation index')
%      xlabel('weeks')

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
    
    case 'make_networkCorr'                           % correlate motor networks for each subj/week
        D       = rs_imana_110118('PP_studySubset'); 
        D       = rs_imana_110118('PP_removePair',D);
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
        save(fullfile(ppDir,'rs_networkCorr.mat'),'-struct','S');          
    case 'make_networkCorr_lesionedIntra'            % correlate motor networks, only lesioned hem in patients
        D       = rs_imana_110118('PP_studySubset'); 
        D       = rs_imana_110118('PP_removePair',D);
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
                Si.mat  = squareform(M);
                S       = addstruct(S,Si);
            end;
        end;
        save(fullfile(ppDir,'rs_networkCorr_lesionedIntra.mat'),'-struct','S');                        
    case 'make_networkCorr_nonlesionedIntra'         % correlate motor networks, only nonlesioned hem in patients
        D       = rs_imana_110118('PP_studySubset'); 
        D       = rs_imana_110118('PP_removePair',D);
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
                Si.mat  = squareform(M);
                S       = addstruct(S,Si);
            end;
        end;
        save(fullfile(ppDir,'rs_networkCorr_nonlesionedIntra.mat'),'-struct','S');                                
    case 'make_all'
    D1 = load(fullfile(ppDir,'rs_networkCorr.mat'));    
    D2 = load(fullfile(ppDir,'rs_networkCorr_lesionedIntra.mat'));
    D3 = load(fullfile(ppDir,'rs_networkCorr_nonlesionedIntra.mat'));
    
    S = D1;
       S.M =[D1.M, D2.M, D3.M];       
       save(fullfile(ppDir,'rs_all.mat'),'-struct','S'); 
    
    case 'plot_networkCorr_all'         % plot correlation across controls and patients
        D   = load(fullfile(ppDir,'rs_networkCorr.mat'));                        
        Dc  = getrow(D,D.control==1);
        Dp  = getrow(D,D.control~=1);
        
        % get avg matrices for each group
        Mc  = fisherinv(nanmean(Dc.M,1));
        Mc  = nonsym_squareform(Mc);
                
        DpW1 = getrow (Dp,Dp.week==1);
        Mp1 = fisherinv(nanmean(DpW1.M,1));
        Mp1 = nonsym_squareform(Mp1);
        DpW4 = getrow (Dp,Dp.week==4);
        Mp4 = fisherinv(nanmean(DpW4.M,1));
        Mp4 = nonsym_squareform(Mp4);
        DpW12 = getrow (Dp,Dp.week==12);
        Mp12 = fisherinv(nanmean(DpW12.M,1));
        Mp12 = nonsym_squareform(Mp12);
        DpW24 = getrow (Dp,Dp.week==24);
        Mp24 = fisherinv(nanmean(DpW24.M,1));
        Mp24 = nonsym_squareform(Mp24);
        DpW52 = getrow (Dp,Dp.week==52);
        Mp52 = fisherinv(nanmean(DpW52.M,1));
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
    
        x1  = fisherinv(nanmean(Dc.mat,1));
        DpW1 = getrow (Dp,Dp.week==1);
        x2 = fisherinv(nanmean(DpW1.mat,1)); 
        DpW4 = getrow (Dp,Dp.week==4);
        x4 = fisherinv(nanmean(DpW4.mat,1));
        DpW12 = getrow (Dp,Dp.week==12);
        x12 = fisherinv(nanmean(DpW12.mat,1));
        DpW24 = getrow (Dp,Dp.week==24);
        x24 = fisherinv(nanmean(DpW24.mat,1));
        DpW52 = getrow (Dp,Dp.week==52);
        x52 = fisherinv(nanmean(DpW52.mat,1));        
        
        Dc = tapply(Dc,{'subj_name'},{'mat','nanmean(x,1)'});
        
        h=make_figure;
        subplot(231);
        title('Week1');
        ylabel('correlation index')
        xlabel('roi pairs')
        traceplot(1:25,fisherinv((Dc.mat)),'errorfcn','stderr');
        hold on
        traceplot(1:25,fisherinv((DpW1.mat)),'errorfcn','stderr','linecolor',[1 0 0],'patchcolor',[1 0 0]);
        ylim([0 1]);        xlim([0 25]);
        subplot(232);
        traceplot(1:25,fisherinv((Dc.mat)),'errorfcn','stderr');
        hold on
        traceplot(1:25,fisherinv((DpW4.mat)),'errorfcn','stderr','linecolor',[1 0 0],'patchcolor',[1 0 0]);
        ylim([0 1]);        xlim([0 25]);
        title('Week4');
        ylabel('correlation index')
        xlabel('roi pairs')
        subplot(233);
        traceplot(1:25,fisherinv((Dc.mat)),'errorfcn','stderr');
        hold on
        traceplot(1:25,fisherinv((DpW12.mat)),'errorfcn','stderr','linecolor',[1 0 0],'patchcolor',[1 0 0]);
        title('Week12');
        ylabel('correlation index')
        xlabel('roi pairs')
        ylim([0 1]);        xlim([0 25]);
        subplot(234);
        traceplot(1:25,fisherinv((Dc.mat)),'errorfcn','stderr');
        hold on
        traceplot(1:25,fisherinv((DpW24.mat)),'errorfcn','stderr','linecolor',[1 0 0],'patchcolor',[1 0 0]);
        title('Week24');
        ylabel('correlation index')
        xlabel('roi pairs')
        ylim([0 1]);        xlim([0 25]);
        subplot(235);
        traceplot(1:25,fisherinv((Dc.mat)),'errorfcn','stderr');
        hold on
        traceplot(1:25,fisherinv((DpW52.mat)),'errorfcn','stderr','linecolor',[1 0 0],'patchcolor',[1 0 0]);
        title('Week52');
        ylim([0 1]);        xlim([0 25]);
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
        
        keyboard;        
    case 'plot_networkCorr_lesioned'    % plot controls and patients (lesioned hemisphere only)
        D   = load(fullfile(ppDir,'rs_networkCorr_lesionedIntra.mat'));        
        
        %take only subgroup
        D   = getrow(D,(D.lesionType~=2)); %no cortical
        D   = getrow(D,(D.lesionType~=4)); %no mixed        
        
        Dc  = getrow(D,D.control==1);
        Dp  = getrow(D,D.control~=1);
        
        % get avg matrices for each group
        Mc  = fisherinv(nanmean(Dc.M,1));
        Mc  = nonsym_squareform(Mc);
                
        DpW1 = getrow (Dp,Dp.week==1);
        Mp1 = fisherinv(nanmean(DpW1.M,1));
        Mp1 = nonsym_squareform(Mp1);
        DpW4 = getrow (Dp,Dp.week==4);
        Mp4 = fisherinv(nanmean(DpW4.M,1));
        Mp4 = nonsym_squareform(Mp4);
        DpW12 = getrow (Dp,Dp.week==12);
        Mp12 = fisherinv(nanmean(DpW12.M,1));
        Mp12 = nonsym_squareform(Mp12);
        DpW24 = getrow (Dp,Dp.week==24);
        Mp24 = fisherinv(nanmean(DpW24.M,1));
        Mp24 = nonsym_squareform(Mp24);
        DpW52 = getrow (Dp,Dp.week==52);
        Mp52 = fisherinv(nanmean(DpW52.M,1));
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
        
         x1  = fisherinv(nanmean(Dc.mat,1));
        DpW1 = getrow (Dp,Dp.week==1);
        x2 = fisherinv(nanmean(DpW1.mat,1)); 
        DpW4 = getrow (Dp,Dp.week==4);
        x4 = fisherinv(nanmean(DpW4.mat,1));
        DpW12 = getrow (Dp,Dp.week==12);
        x12 = fisherinv(nanmean(DpW12.mat,1));
        DpW24 = getrow (Dp,Dp.week==24);
        x24 = fisherinv(nanmean(DpW24.mat,1));
        DpW52 = getrow (Dp,Dp.week==52);
        x52 = fisherinv(nanmean(DpW52.mat,1));        
        
        Dc = tapply(Dc,{'subj_name'},{'mat','nanmean(x,1)'});
        
        h=make_figure;
        subplot(231);
        title('Week1');
        ylabel('correlation index')
        xlabel('roi pairs')
        traceplot(1:10,fisherinv((Dc.mat)),'errorfcn','stderr');
        hold on
        traceplot(1:10,fisherinv((DpW1.mat)),'errorfcn','stderr','linecolor',[1 0 0],'patchcolor',[1 0 0]);
        ylim([0 1]);        xlim([0 10]);
        subplot(232);
        traceplot(1:10,fisherinv((Dc.mat)),'errorfcn','stderr');
        hold on
        traceplot(1:10,fisherinv((DpW4.mat)),'errorfcn','stderr','linecolor',[1 0 0],'patchcolor',[1 0 0]);
        ylim([0 1]);        xlim([0 10]);
        title('Week4');
        ylabel('correlation index')
        xlabel('roi pairs')
        subplot(233);
        traceplot(1:10,fisherinv((Dc.mat)),'errorfcn','stderr');
        hold on
        traceplot(1:10,fisherinv((DpW12.mat)),'errorfcn','stderr','linecolor',[1 0 0],'patchcolor',[1 0 0]);
        title('Week12');
        ylabel('correlation index')
        xlabel('roi pairs')
        ylim([0 1]);        xlim([0 10]);
        subplot(234);
        traceplot(1:10,fisherinv((Dc.mat)),'errorfcn','stderr');
        hold on
        traceplot(1:10,fisherinv((DpW24.mat)),'errorfcn','stderr','linecolor',[1 0 0],'patchcolor',[1 0 0]);
        title('Week24');
        ylabel('correlation index')
        xlabel('roi pairs')
        ylim([0 1]);        xlim([0 10]);
        subplot(235);
        traceplot(1:10,fisherinv((Dc.mat)),'errorfcn','stderr');
        hold on
        traceplot(1:10,fisherinv((DpW52.mat)),'errorfcn','stderr','linecolor',[1 0 0],'patchcolor',[1 0 0]);
        title('Week52');
        ylim([0 1]);        xlim([0 10]);
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
        
        keyboard;        
    case 'plot_networkCorr_nonlesioned' % plot controls and patients (non lesioned hemisphere only)
        D   = load(fullfile(ppDir,'rs_networkCorr_nonlesionedIntra.mat'));        
        
        %take only subgroup
        D   = getrow(D,(D.lesionType~=2)); %no cortical
        D   = getrow(D,(D.lesionType~=4)); %no mixed

        Dc  = getrow(D,D.control==1);
        Dp  = getrow(D,D.control~=1);
        
        % get avg matrices for each group
        Mc  = fisherinv(nanmean(Dc.M,1));
        Mc  = nonsym_squareform(Mc);
        
        
        DpW1 = getrow (Dp,Dp.week==1);
        Mp1 = fisherinv(nanmean(DpW1.M,1));
        Mp1 = nonsym_squareform(Mp1);
        DpW4 = getrow (Dp,Dp.week==4);
        Mp4 = fisherinv(nanmean(DpW4.M,1));
        Mp4 = nonsym_squareform(Mp4);
        DpW12 = getrow (Dp,Dp.week==12);
        Mp12 = fisherinv(nanmean(DpW12.M,1));
        Mp12 = nonsym_squareform(Mp12);
        DpW24 = getrow (Dp,Dp.week==24);
        Mp24 = fisherinv(nanmean(DpW24.M,1));
        Mp24 = nonsym_squareform(Mp24);
        DpW52 = getrow (Dp,Dp.week==52);
        Mp52 = fisherinv(nanmean(DpW52.M,1));
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
    
 x1  = fisherinv(nanmean(Dc.mat,1));
        DpW1 = getrow (Dp,Dp.week==1);
        x2 = fisherinv(nanmean(DpW1.mat,1)); 
        DpW4 = getrow (Dp,Dp.week==4);
        x4 = fisherinv(nanmean(DpW4.mat,1));
        DpW12 = getrow (Dp,Dp.week==12);
        x12 = fisherinv(nanmean(DpW12.mat,1));
        DpW24 = getrow (Dp,Dp.week==24);
        x24 = fisherinv(nanmean(DpW24.mat,1));
        DpW52 = getrow (Dp,Dp.week==52);
        x52 = fisherinv(nanmean(DpW52.mat,1));        
        
        Dc = tapply(Dc,{'subj_name'},{'mat','nanmean(x,1)'});
        
        h=make_figure;
        subplot(231);
        title('Week1');
        ylabel('correlation index')
        xlabel('roi pairs')
        traceplot(1:10,fisherinv((Dc.mat)),'errorfcn','stderr');
        hold on
        traceplot(1:10,fisherinv((DpW1.mat)),'errorfcn','stderr','linecolor',[1 0 0],'patchcolor',[1 0 0]);
        ylim([0 1]);        xlim([0 10]);
        subplot(232);
        traceplot(1:10,fisherinv((Dc.mat)),'errorfcn','stderr');
        hold on
        traceplot(1:10,fisherinv((DpW4.mat)),'errorfcn','stderr','linecolor',[1 0 0],'patchcolor',[1 0 0]);
        ylim([0 1]);        xlim([0 10]);
        title('Week4');
        ylabel('correlation index')
        xlabel('roi pairs')
        subplot(233);
        traceplot(1:10,fisherinv((Dc.mat)),'errorfcn','stderr');
        hold on
        traceplot(1:10,fisherinv((DpW12.mat)),'errorfcn','stderr','linecolor',[1 0 0],'patchcolor',[1 0 0]);
        title('Week12');
        ylabel('correlation index')
        xlabel('roi pairs')
        ylim([0 1]);        xlim([0 10]);
        subplot(234);
        traceplot(1:10,fisherinv((Dc.mat)),'errorfcn','stderr');
        hold on
        traceplot(1:10,fisherinv((DpW24.mat)),'errorfcn','stderr','linecolor',[1 0 0],'patchcolor',[1 0 0]);
        title('Week24');
        ylabel('correlation index')
        xlabel('roi pairs')
        ylim([0 1]);        xlim([0 10]);
        subplot(235);
        traceplot(1:10,fisherinv((Dc.mat)),'errorfcn','stderr');
        hold on
        traceplot(1:10,fisherinv((DpW52.mat)),'errorfcn','stderr','linecolor',[1 0 0],'patchcolor',[1 0 0]);
        title('Week52');
        ylim([0 1]);        xlim([0 10]);
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
        
        keyboard;                     
    case 'plot_all_conn'
     D = load(fullfile(ppDir,'nc_getThemall.mat'));                       
        Dc  = getrow(D,D.control==1);
        Dp  = getrow(D,D.control~=1);
        
        % get avg matrices for each group
        Mc  = fisherinv(nanmean(Dc.M,1));
        Mc  = squareform(Mc);
        
                
        DpW1 = getrow (Dp,Dp.week==1);
        Mp1 = fisherinv(nanmean(DpW1.M,1));
        Mp1 = squareform(Mp1);
        DpW4 = getrow (Dp,Dp.week==4);
        Mp4 = fisherinv(nanmean(DpW4.M,1));
        Mp4 = squareform(Mp4);
        DpW12 = getrow (Dp,Dp.week==12);
        Mp12 = fisherinv(nanmean(DpW12.M,1));
        Mp12 = squareform(Mp12);
        DpW24 = getrow (Dp,Dp.week==24);
        Mp24 = fisherinv(nanmean(DpW24.M,1));
        Mp24 = squareform(Mp24);
        DpW52 = getrow (Dp,Dp.week==52);
        Mp52 = fisherinv(nanmean(DpW52.M,1));
        Mp52 = squareform(Mp52);
        
        
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
    
        x1  = fisherinv(nanmean(Dc.M,1));
        DpW1 = getrow (Dp,Dp.week==1);
        x2 = fisherinv(nanmean(DpW1.M,1)); 
        DpW4 = getrow (Dp,Dp.week==4);
        x4 = fisherinv(nanmean(DpW4.M,1));
        DpW12 = getrow (Dp,Dp.week==12);
        x12 = fisherinv(nanmean(DpW12.M,1));
        DpW24 = getrow (Dp,Dp.week==24);
        x24 = fisherinv(nanmean(DpW24.M,1));
        DpW52 = getrow (Dp,Dp.week==52);
        x52 = fisherinv(nanmean(DpW52.M,1));        
        
        Dc = tapply(Dc,{'subj_name'},{'M','nanmean(x,1)'});
        
        h=make_figure;
        subplot(231);
        title('Week1');
        ylabel('correlation index')
        xlabel('roi pairs')
        traceplot(1:45,fisherinv((Dc.M)),'errorfcn','stderr');
        hold on
        traceplot(1:45,fisherinv((DpW1.M)),'errorfcn','stderr','linecolor',[1 0 0],'patchcolor',[1 0 0]);
        ylim([0 1]);        xlim([0 45]);
        subplot(232);
        traceplot(1:45,fisherinv((Dc.M)),'errorfcn','stderr');
        hold on
        traceplot(1:45,fisherinv((DpW4.M)),'errorfcn','stderr','linecolor',[1 0 0],'patchcolor',[1 0 0]);
        ylim([0 1]);        xlim([0 45]);
        title('Week4');
        ylabel('correlation index')
        xlabel('roi pairs')
        subplot(233);
        traceplot(1:45,fisherinv((Dc.M)),'errorfcn','stderr');
        hold on
        traceplot(1:45,fisherinv((DpW12.M)),'errorfcn','stderr','linecolor',[1 0 0],'patchcolor',[1 0 0]);
        title('Week12');
        ylabel('correlation index')
        xlabel('roi pairs')
        ylim([0 1]);        xlim([0 45]);
        subplot(234);
        traceplot(1:45,fisherinv((Dc.M)),'errorfcn','stderr');
        hold on
        traceplot(1:45,fisherinv((DpW24.M)),'errorfcn','stderr','linecolor',[1 0 0],'patchcolor',[1 0 0]);
        title('Week24');
        ylabel('correlation index')
        xlabel('roi pairs')
        ylim([0 1]);        xlim([0 45]);
        subplot(235);
        traceplot(1:45,fisherinv((Dc.M)),'errorfcn','stderr');
        hold on
        traceplot(1:45,fisherinv((DpW52.M)),'errorfcn','stderr','linecolor',[1 0 0],'patchcolor',[1 0 0]);
        title('Week52');
        ylim([0 1]);        xlim([0 45]);
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
        
    case 'MainFigure_Inter'
        D   = load(fullfile(ppDir,'nc_interhem.mat'));                
        Dc  = getrow(D,D.control==1);
        Dp  = getrow(D,D.control==0);
        
        Dc  = tapply(Dc,{'subj_name'},{'M','mean(x,1)'});
        Dp  = tapply(Dp,{'week'},{'M','mean(x,1)'});
        
        C = [];
        P = [];
        for c=unique(Dc.subj_name)'
            isC     = getrow(Dc,ismember(Dc.subj_name,c));
            notC    = getrow(Dc,~ismember(Dc.subj_name,c));
            
            % add control correlations with controls
            Ci.subj_name    = c;
            Ci.week         = 0;
            Ci.fz_r         = fisherz(corr(isC.M',mean(notC.M,1)'));
            C               = addstruct(C,Ci);
            
            % add control correlations with patients
            Pi.subj_name    = repmat(c,5,1);
            Pi.week         = [1;4;12;24;52];
            Pi.fz_r         = fisherz(corr(isC.M',Dp.M'))';
            P               = addstruct(P,Pi);
        end;
        
        %stats
        keyboard;
        
        
        % make plots
        CAT.linecolor   = 'r';
        CAT.errorcolor  = 'r';
        CAT.markerfill  = 'r';
        CAT.markersize  = 8;
        
        lineplot(P.week,P.fz_r,'plotfcn','fisherinv(mean(x))','CAT',CAT);
        m = mean(C.fz_r); se = stderr(C.fz_r);
        drawline(fisherinv(m),'dir','horz','color','b');
        drawline(fisherinv(m+se),'dir','horz','color','b','linestyle','--');
        drawline(fisherinv(m-se),'dir','horz','color','b','linestyle','--');
        title('interhem')
    case 'MainFigure_Intra_Lesioned'
        D   = load(fullfile(ppDir,'nc_intrahem_lesioned.mat'));                
        Dc  = getrow(D,D.control==1);
        Dp  = getrow(D,D.control==0);
        
        Dc  = tapply(Dc,{'subj_name'},{'M','mean(x,1)'});
        Dp  = tapply(Dp,{'week'},{'M','mean(x,1)'});
        
        C = [];
        P = [];
        for c=unique(Dc.subj_name)'
            isC     = getrow(Dc,ismember(Dc.subj_name,c));
            notC    = getrow(Dc,~ismember(Dc.subj_name,c));
            
            % add control correlations with controls
            Ci.subj_name    = c;
            Ci.week         = 0;
            Ci.fz_r         = fisherz(corr(isC.M',mean(notC.M,1)'));
            C               = addstruct(C,Ci);
            
            % add control correlations with patients
            Pi.subj_name    = repmat(c,5,1);
            Pi.week         = [1;4;12;24;52];
            Pi.fz_r         = fisherz(corr(isC.M',Dp.M'))';
            P               = addstruct(P,Pi);
        end;
        
        % make plots
        CAT.linecolor   = 'r';
        CAT.errorcolor  = 'r';
        CAT.markerfill  = 'r';
        CAT.markersize  = 8;
        
        lineplot(P.week,P.fz_r,'plotfcn','fisherinv(mean(x))','CAT',CAT);
        m = mean(C.fz_r); se = stderr(C.fz_r);
        drawline(fisherinv(m),'dir','horz','color','b');
        drawline(fisherinv(m+se),'dir','horz','color','b','linestyle','--');
        drawline(fisherinv(m-se),'dir','horz','color','b','linestyle','--');
        title('intralesioned')       
    case 'MainFigure_Intra_NonLesioned'
        D   = load(fullfile(ppDir,'nc_intrahem_nonlesioned.mat'));                
        Dc  = getrow(D,D.control==1);
        Dp  = getrow(D,D.control==0);
        
        Dc  = tapply(Dc,{'subj_name'},{'M','mean(x,1)'});
        Dp  = tapply(Dp,{'week'},{'M','mean(x,1)'});
        
        C = [];
        P = [];
        for c=unique(Dc.subj_name)'
            isC     = getrow(Dc,ismember(Dc.subj_name,c));
            notC    = getrow(Dc,~ismember(Dc.subj_name,c));
            
            % add control correlations with controls
            Ci.subj_name    = c;
            Ci.week         = 0;
            Ci.fz_r         = fisherz(corr(isC.M',mean(notC.M,1)'));
            C               = addstruct(C,Ci);
            
            % add control correlations with patients
            Pi.subj_name    = repmat(c,5,1);
            Pi.week         = [1;4;12;24;52];
            Pi.fz_r         = fisherz(corr(isC.M',Dp.M'))';
            P               = addstruct(P,Pi);
        end;
        
        % make plots
        CAT.linecolor   = 'r';
        CAT.errorcolor  = 'r';
        CAT.markerfill  = 'r';
        CAT.markersize  = 8;
        
        lineplot(P.week,P.fz_r,'plotfcn','fisherinv(mean(x))','CAT',CAT);
        m = mean(C.fz_r); se = stderr(C.fz_r);
        drawline(fisherinv(m),'dir','horz','color','b');
        drawline(fisherinv(m+se),'dir','horz','color','b','linestyle','--');
        drawline(fisherinv(m-se),'dir','horz','color','b','linestyle','--');
        title('intranonlesioned')               
    case 'MainFigure_all_conn'
        D = load(fullfile(ppDir,'nc_getThemall.mat'));
        
        Dc  = getrow(D,D.control==1);
        Dp  = getrow(D,D.control==0);
        
        Dc  = tapply(Dc,{'subj_name'},{'M','mean(x,1)'});
        Dp  = tapply(Dp,{'week'},{'M','mean(x,1)'});
        
        C = [];
        P = [];
        for c=unique(Dc.subj_name)'
            isC     = getrow(Dc,ismember(Dc.subj_name,c));
            notC    = getrow(Dc,~ismember(Dc.subj_name,c));
            
            % add control correlations with controls
            Ci.subj_name    = c;
            Ci.week         = 0;
            Ci.fz_r         = fisherz(corr(isC.M',mean(notC.M,1)'));
            C               = addstruct(C,Ci);
            
            % add control correlations with patients
            Pi.subj_name    = repmat(c,5,1);
            Pi.week         = [1;4;12;24;52];
            Pi.fz_r         = fisherz(corr(isC.M',Dp.M'))';
            P               = addstruct(P,Pi);
        end;
        
        %stats
        keyboard;
        [W,V] = pivottable(P.subj_name,P.week,P.fz_r,'mean')
        %take W for statistics in spss
        
        % make plots
        CAT.linecolor   = 'r';
        CAT.errorcolor  = 'r';
        CAT.markerfill  = 'r';
        CAT.markersize  = 8;
        
        lineplot(P.week,P.fz_r,'plotfcn','fisherinv(mean(x))','CAT',CAT);
        m = mean(C.fz_r); se = stderr(C.fz_r);
        drawline(fisherinv(m),'dir','horz','color','b');
        drawline(fisherinv(m+se),'dir','horz','color','b','linestyle','--');
        drawline(fisherinv(m-se),'dir','horz','color','b','linestyle','--');
        title('interhem')
        
    case 'MainFigure_All'
        subplot(221);
        rs_imana_110118('MainFigure_Inter');
        ylim([0.5 1.1]);
        ylabel('Fisher Z')
        xlabel('weeks')
        subplot(222);
        rs_imana_110118('MainFigure_Intra_Lesioned');
        ylim([0.5 1.1]);
         ylabel('Fisher Z')
        xlabel('weeks')
        subplot(223);
        rs_imana_110118('MainFigure_Intra_NonLesioned');
        ylim([0.5 1.1]); 
        ylabel('Fisher Z')
        xlabel('weeks')        
        subplot(224);
        rs_imana_110118('MainFigure_all_conn');
        ylim([0.5 1.1]); 
        ylabel('Fisher Z')
        xlabel('weeks')   
        
    case 'MainFigure_Inter_EUC'
        D   = load(fullfile(ppDir,'nc_interhem.mat'));                
        Dc  = getrow(D,D.control==1);
        Dp  = getrow(D,D.control==0);
        
        Dc  = tapply(Dc,{'subj_name'},{'M','mean(x,1)'});
        Dp  = tapply(Dp,{'week'},{'M','mean(x,1)'});
        
        C = [];
        P = [];
        leaveoutCAvg = nan(length(Dc.subj_name),25);
        

        idx = 1;
        for c=unique(Dc.subj_name)'
            isC     = getrow(Dc,ismember(Dc.subj_name,c));
            notC    = getrow(Dc,~ismember(Dc.subj_name,c));
            
            % add control correlations with controls
            Ci.subj_name    = c;
            Ci.week         = 0;
            
            % calculate leave-one-out average of all controls
            leaveoutCAvg(idx,:)  = nanmean(notC.M,1);
            Ci.euc              = pdist2(isC.M,leaveoutCAvg(idx,:),'euclidean');
            
            idx             = idx + 1;
            C               = addstruct(C,Ci);
            
            % add control correlations with patients
            for j = 1:length(Dp.week)
            Pi.subj_name    = c;
            Pi.week         = Dp.week(j,:);
            Pi.euc          = pdist2(isC.M,Dp.M(j,:),'euclidean');
            P               = addstruct(P,Pi);
            end
            
        end;  
      keyboard;
      dsave(fullfile(statsDir,'MainFigure_Inter_EUC_P.dat'),P); 
      dsave(fullfile(statsDir,'MainFigure_Inter_EUC_C.dat'),C); 
      
        % make plots
        CAT.linecolor   = 'r';
        CAT.errorcolor  = 'r';
        CAT.markerfill  = 'r';
        CAT.markersize  = 8;
        
        lineplot(P.week,P.euc,'plotfcn','(mean(x))','CAT',CAT);
        m = mean(C.euc); se = stderr(C.euc);
        drawline(m,'dir','horz','color','b');
        drawline(m+se,'dir','horz','color','b','linestyle','--');
        drawline(m-se,'dir','horz','color','b','linestyle','--');
        title('interhem')
        keyboard;
    case 'MainFigure_Intra_Lesioned_EUC'
        D   = load(fullfile(ppDir,'nc_intrahem_lesioned.mat'));                                
        Dc  = getrow(D,D.control==1);
        Dp  = getrow(D,D.control==0);
        
        Dc  = tapply(Dc,{'subj_name'},{'M','mean(x,1)'});
        Dp  = tapply(Dp,{'week'},{'M','mean(x,1)'});
        
        C = [];
        P = [];
        leaveoutCAvg = nan(length(Dc.subj_name),10);
        

        idx = 1;
        for c=unique(Dc.subj_name)'
            isC     = getrow(Dc,ismember(Dc.subj_name,c));
            notC    = getrow(Dc,~ismember(Dc.subj_name,c));
            
            % add control correlations with controls
            Ci.subj_name    = c;
            Ci.week         = 0;
            
            % calculate leave-one-out average of all controls
            leaveoutCAvg(idx,:)  = nanmean(notC.M,1);
            Ci.euc              = pdist2(isC.M,leaveoutCAvg(idx,:),'euclidean');
            
            idx             = idx + 1;
            C               = addstruct(C,Ci);
            
            % add control correlations with patients
            for j = 1:length(Dp.week)
            Pi.subj_name    = c;
            Pi.week         = Dp.week(j,:);
            Pi.euc          = pdist2(isC.M,Dp.M(j,:),'euclidean');
            P               = addstruct(P,Pi);
            end
            
        end;  
      keyboard;
      dsave(fullfile(statsDir,'MainFigure_Intra_Lesioned_EUC_P.dat'),P); 
      dsave(fullfile(statsDir,'MainFigure_Intra_Lesioned_EUC_C.dat'),C); 
      
      
        % make plots
        CAT.linecolor   = 'r';
        CAT.errorcolor  = 'r';
        CAT.markerfill  = 'r';
        CAT.markersize  = 8;
        
        lineplot(P.week,P.euc,'plotfcn','(mean(x))','CAT',CAT);
        m = mean(C.euc); se = stderr(C.euc);
        drawline(m,'dir','horz','color','b');
        drawline(m+se,'dir','horz','color','b','linestyle','--');
        drawline(m-se,'dir','horz','color','b','linestyle','--');
        title('intrahem lesioned')      
    case 'MainFigure_Intra_NonLesioned_EUC'
        D   = load(fullfile(ppDir,'nc_intrahem_nonlesioned.mat'));                                               
        Dc  = getrow(D,D.control==1);
        Dp  = getrow(D,D.control==0);
        
        Dc  = tapply(Dc,{'subj_name'},{'M','mean(x,1)'});
        Dp  = tapply(Dp,{'week'},{'M','mean(x,1)'});
        
        C = [];
        P = [];
        leaveoutCAvg = nan(length(Dc.subj_name),10);
        

        idx = 1;
        for c=unique(Dc.subj_name)'
            isC     = getrow(Dc,ismember(Dc.subj_name,c));
            notC    = getrow(Dc,~ismember(Dc.subj_name,c));
            
            % add control correlations with controls
            Ci.subj_name    = c;
            Ci.week         = 0;
            
            % calculate leave-one-out average of all controls
            leaveoutCAvg(idx,:)  = nanmean(notC.M,1);
            Ci.euc              = pdist2(isC.M,leaveoutCAvg(idx,:),'euclidean');
            
            idx             = idx + 1;
            C               = addstruct(C,Ci);
            
            % add control correlations with patients
            for j = 1:length(Dp.week)
            Pi.subj_name    = c;
            Pi.week         = Dp.week(j,:);
            Pi.euc          = pdist2(isC.M,Dp.M(j,:),'euclidean');
            P               = addstruct(P,Pi);
            end
            
        end;  
      keyboard;
      
      dsave(fullfile(statsDir,'MainFigure_Intra_NonLesioned_EUC_P.dat'),P); 
      dsave(fullfile(statsDir,'MainFigure_Intra_NonLesioned_EUC_C.dat'),C);   
      
      % make plots
        CAT.linecolor   = 'r';
        CAT.errorcolor  = 'r';
        CAT.markerfill  = 'r';
        CAT.markersize  = 8;
        
        lineplot(P.week,P.euc,'plotfcn','(mean(x))','CAT',CAT);
        m = mean(C.euc); se = stderr(C.euc);
        drawline(m,'dir','horz','color','b');
        drawline(m+se,'dir','horz','color','b','linestyle','--');
        drawline(m-se,'dir','horz','color','b','linestyle','--');
        title('intranonlesioned')  
    case 'MainFigure_all_conn_EUC'
       D = load(fullfile(ppDir,'nc_getThemall.mat'));
        Dc  = getrow(D,D.control==1);
        Dp  = getrow(D,D.control==0);
        
        Dc  = tapply(Dc,{'subj_name'},{'M','mean(x,1)'});
        Dp  = tapply(Dp,{'week'},{'M','mean(x,1)'});
        
        C = [];
        P = [];
        leaveoutCAvg = nan(length(Dc.subj_name),45);
        

        idx = 1;
        for c=unique(Dc.subj_name)'
            isC     = getrow(Dc,ismember(Dc.subj_name,c));
            notC    = getrow(Dc,~ismember(Dc.subj_name,c));
            
            % add control correlations with controls
            Ci.subj_name    = c;
            Ci.week         = 0;
            
            % calculate leave-one-out average of all controls
            leaveoutCAvg(idx,:)  = nanmean(notC.M,1);
            Ci.euc              = pdist2(isC.M,leaveoutCAvg(idx,:),'euclidean');
            
            idx             = idx + 1;
            C               = addstruct(C,Ci);
            
            % add control correlations with patients
            for j = 1:length(Dp.week)
            Pi.subj_name    = c;
            Pi.week         = Dp.week(j,:);
            Pi.euc          = pdist2(isC.M,Dp.M(j,:),'euclidean');
            P               = addstruct(P,Pi);
            end
            
        end;  
      keyboard;
      dsave(fullfile(statsDir,'MainFigure_all_EUC_P.dat'),P); 
      dsave(fullfile(statsDir,'MainFigure_all_EUC_C.dat'),C); 
      
        % make plots
        CAT.linecolor   = 'r';
        CAT.errorcolor  = 'r';
        CAT.markerfill  = 'r';
        CAT.markersize  = 8;
        
        lineplot(P.week,P.euc,'plotfcn','(mean(x))','CAT',CAT);
        m = mean(C.euc); se = stderr(C.euc);
        drawline(m,'dir','horz','color','b');
        drawline(m+se,'dir','horz','color','b','linestyle','--');
        drawline(m-se,'dir','horz','color','b','linestyle','--');
        title('all connections')  
    case 'MainFigure_EUC'
        subplot(221);
        rs_imana_110118('MainFigure_Inter_EUC');
        ylim([0.5 3]);
        ylabel('Euclidian distance')
        xlabel('weeks')
        subplot(222);
        rs_imana_110118('MainFigure_Intra_Lesioned_EUC');
        ylim([0.5 3]);
         ylabel('Euclidian distance')
        xlabel('weeks')
        subplot(223);
        rs_imana_110118('MainFigure_Intra_NonLesioned_EUC');
        ylim([0.5 3]); 
        ylabel('Euclidian distance')
        xlabel('weeks')        
        subplot(224);
        rs_imana_110118('MainFigure_all_conn_EUC');
        ylim([0.5 3]); 
        ylabel('Euclidian distance')
        xlabel('weeks')        
        keyboard;

        
        
    case 'MainFigure_Inter_EUC_patients'
        D   = load(fullfile(ppDir,'nc_interhem.mat'));                                             
        Dc  = getrow(D,D.control==1);
        Dp  = getrow(D,D.control==0);
        
        Dc  = tapply(Dc,{'subj_name'},{'M','mean(x,1)'});
        Dp  = tapply(Dp,{'subj_name','week'},{'M','mean(x,1)'});
        
        C = [];
        P = [];
        leaveoutCAvg = nan(length(Dc.subj_name),25);
        Cavg         = nanmean(Dc.M);
        
        idx = 1;
        for c=unique(Dc.subj_name)'
            isC     = getrow(Dc,ismember(Dc.subj_name,c));
            notC    = getrow(Dc,~ismember(Dc.subj_name,c));
            
            % add control correlations with controls
            Ci.subj_name    = c;
            Ci.week         = 0;
            
            % calculate leave-one-out average of all controls
            leaveoutCAvg(idx,:)  = nanmean(notC.M,1);
            Ci.euc          = pdist2(isC.M,leaveoutCAvg(idx,:),'euclidean');
            
            idx             = idx + 1;
            C               = addstruct(C,Ci);
        end;
       
        
        % calculate euclidean distance for patients to avg control
        for p=unique(Dp.subj_name)'
            isP     = getrow(Dp,ismember(Dp.subj_name,p));
            for j = 1:length(isP.week)
                Pi.subj_name  = p;
                Pi.week       = isP.week(j,:);          
                Pi.euc  = pdist2(isP.M(j,:),Cavg,'euclidean');
                P       = addstruct(P,Pi);
            end 
            
        end
      dsave(fullfile(statsDir,'MainFigure_Inter_EUC_patients_P.dat'),P); 
      dsave(fullfile(statsDir,'MainFigure_Inter_EUC_patients_C.dat'),C); 
      
        % make plots
        CAT.linecolor   = 'r';
        CAT.errorcolor  = 'r';
        CAT.markerfill  = 'r';
        CAT.markersize  = 8;
        
        lineplot(P.week,P.euc,'plotfcn','(mean(x))','CAT',CAT);
        m = mean(C.euc); se = stderr(C.euc);
        drawline(m,'dir','horz','color','b');
        drawline(m+se,'dir','horz','color','b','linestyle','--');
        drawline(m-se,'dir','horz','color','b','linestyle','--');
        title('interhem')
        keyboard;
    case 'MainFigure_Intra_Lesioned_EUC_patients'
        D   = load(fullfile(ppDir,'nc_intrahem_lesioned.mat'));    
        Dc  = getrow(D,D.control==1);
        Dp  = getrow(D,D.control==0);
        
        Dc  = tapply(Dc,{'subj_name'},{'M','mean(x,1)'});
        Dp  = tapply(Dp,{'subj_name','week'},{'M','mean(x,1)'});
        
        C = [];
        P = [];
        
        leaveoutCAvg = nan(length(Dc.subj_name),10);
        Cavg         = nanmean(Dc.M);
        
        idx = 1;
        for c=unique(Dc.subj_name)'
            isC     = getrow(Dc,ismember(Dc.subj_name,c));
            notC    = getrow(Dc,~ismember(Dc.subj_name,c));
            
            % add control correlations with controls
            Ci.subj_name    = c;
            Ci.week         = 0;
            
            % calculate leave-one-out average of all controls
            leaveoutCAvg(idx,:)  = nanmean(notC.M,1);
            Ci.euc          = pdist2(isC.M,leaveoutCAvg(idx,:),'euclidean');
            
            idx             = idx + 1;
            C               = addstruct(C,Ci);
        end;
        
         % calculate euclidean distance for patients to avg control
        for p=unique(Dp.subj_name)'
            isP     = getrow(Dp,ismember(Dp.subj_name,p));
            for j = 1:length(isP.week)
                Pi.subj_name  = p;
                Pi.week       = isP.week(j,:);          
                Pi.euc  = pdist2(isP.M(j,:),Cavg,'euclidean');
                P       = addstruct(P,Pi);
            end 
            
        end
        
      dsave(fullfile(statsDir,'MainFigure_Intra_Lesioned_EUC_patients_P.dat'),P); 
      dsave(fullfile(statsDir,'MainFigure_Intra_Lesioned_EUC_patients_C.dat'),C); 
        
        % make plots
        CAT.linecolor   = 'r';
        CAT.errorcolor  = 'r';
        CAT.markerfill  = 'r';
        CAT.markersize  = 8;
        
        lineplot(P.week,P.euc,'plotfcn','(mean(x))','CAT',CAT);
        m = mean(C.euc); se = stderr(C.euc);
        drawline(m,'dir','horz','color','b');
        drawline(m+se,'dir','horz','color','b','linestyle','--');
        drawline(m-se,'dir','horz','color','b','linestyle','--');
        title('intra lesioned')
    case 'MainFigure_Intra_NonLesioned_EUC_patients'  
        D   = load(fullfile(ppDir,'nc_intrahem_nonlesioned.mat'));                                            
        Dc  = getrow(D,D.control==1);
        Dp  = getrow(D,D.control==0);
        
        Dc  = tapply(Dc,{'subj_name'},{'M','mean(x,1)'});
        Dp  = tapply(Dp,{'subj_name','week'},{'M','mean(x,1)'});
        
        C = [];
        P = [];
        leaveoutCAvg = nan(length(Dc.subj_name),10);
        Cavg         = nanmean(Dc.M);
        
        idx = 1;
        for c=unique(Dc.subj_name)'
            isC     = getrow(Dc,ismember(Dc.subj_name,c));
            notC    = getrow(Dc,~ismember(Dc.subj_name,c));
            
            % add control correlations with controls
            Ci.subj_name    = c;
            Ci.week         = 0;
            
            % calculate leave-one-out average of all controls
            leaveoutCAvg(idx,:)  = nanmean(notC.M,1);
            Ci.euc          = pdist2(isC.M,leaveoutCAvg(idx,:),'euclidean');
            
            idx             = idx + 1;
            C               = addstruct(C,Ci);
        end;
        
         % calculate euclidean distance for patients to avg control
        for p=unique(Dp.subj_name)'
            isP     = getrow(Dp,ismember(Dp.subj_name,p));
            for j = 1:length(isP.week)
                Pi.subj_name  = p;
                Pi.week       = isP.week(j,:);          
                Pi.euc  = pdist2(isP.M(j,:),Cavg,'euclidean');
                P       = addstruct(P,Pi);
            end 
            
        end
        
      dsave(fullfile(statsDir,'MainFigure_Intra_NonLesioned_EUC_patients_P.dat'),P); 
      dsave(fullfile(statsDir,'MainFigure_Intra_NonLesioned_EUC_patients_C.dat'),C); 
        
        % make plots
        CAT.linecolor   = 'r';
        CAT.errorcolor  = 'r';
        CAT.markerfill  = 'r';
        CAT.markersize  = 8;
        
        lineplot(P.week,P.euc,'plotfcn','(mean(x))','CAT',CAT);
        m = mean(C.euc); se = stderr(C.euc);
        drawline(m,'dir','horz','color','b');
        drawline(m+se,'dir','horz','color','b','linestyle','--');
        drawline(m-se,'dir','horz','color','b','linestyle','--');
        title('intra nonlesioned')
    case 'MainFigure_all_conn_EUC_patients' 
        D = load(fullfile(ppDir,'nc_getThemall.mat'));                                               
        Dc  = getrow(D,D.control==1);
        Dp  = getrow(D,D.control==0);
        
        Dc  = tapply(Dc,{'subj_name'},{'M','mean(x,1)'});
        Dp  = tapply(Dp,{'subj_name','week'},{'M','mean(x,1)'});
        
        C = [];
        P = [];
        leaveoutCAvg = nan(length(Dc.subj_name),45);
        Cavg         = nanmean(Dc.M);
        
        idx = 1;
        for c=unique(Dc.subj_name)'
            isC     = getrow(Dc,ismember(Dc.subj_name,c));
            notC    = getrow(Dc,~ismember(Dc.subj_name,c));
            
            % add control correlations with controls
            Ci.subj_name    = c;
            Ci.week         = 0;
            
            % calculate leave-one-out average of all controls
            leaveoutCAvg(idx,:)  = nanmean(notC.M,1);
            Ci.euc          = pdist2(isC.M,leaveoutCAvg(idx,:),'euclidean');
            
            idx             = idx + 1;
            C               = addstruct(C,Ci);
        end;
        
         % calculate euclidean distance for patients to avg control
        for p=unique(Dp.subj_name)'
            isP     = getrow(Dp,ismember(Dp.subj_name,p));
            for j = 1:length(isP.week)
                Pi.subj_name  = p;
                Pi.week       = isP.week(j,:);          
                Pi.euc  = pdist2(isP.M(j,:),Cavg,'euclidean');
                P       = addstruct(P,Pi);
            end 
            
        end
        
      dsave(fullfile(statsDir,'MainFigure_all_conn_EUC_patients_P.dat'),P); 
      dsave(fullfile(statsDir,'MainFigure_all_conn_EUC_patients_C.dat'),C); 
        
        % make plots
        CAT.linecolor   = 'r';
        CAT.errorcolor  = 'r';
        CAT.markerfill  = 'r';
        CAT.markersize  = 8;
        
        lineplot(P.week,P.euc,'plotfcn','(mean(x))','CAT',CAT);
        m = mean(C.euc); se = stderr(C.euc);
        drawline(m,'dir','horz','color','b');
        drawline(m+se,'dir','horz','color','b','linestyle','--');
        drawline(m-se,'dir','horz','color','b','linestyle','--');
        title('all connections')    
    case 'MainFigure_EUC_patients' 
        subplot(221);
        rs_imana_110118('MainFigure_Inter_EUC_patients');
        ylim([0.5 4]);
        ylabel('Euclidian distance')
        xlabel('weeks')
        subplot(222);
        rs_imana_110118('MainFigure_Intra_Lesioned_EUC_patients');
        ylim([0.5 4]);
         ylabel('Euclidian distance')
        xlabel('weeks')
        subplot(223);
        rs_imana_110118('MainFigure_Intra_NonLesioned_EUC_patients');
        ylim([0.5 4]); 
        ylabel('Euclidian distance')
        xlabel('weeks')        
        subplot(224);
        rs_imana_110118('MainFigure_all_conn_EUC_patients');
        ylim([0.5 4]); 
        ylabel('Euclidian distance')
        xlabel('weeks')        
        keyboard;
        
    case 'residual_analysis'
        % build mean vector for controls 
        D   = load(fullfile(ppDir,'nc_interhem.mat'));                
        Dc  = getrow(D,D.control==1);
        Dp  = getrow(D,D.control==0);
        
        
                c = mean(Dc.M,1);
                Dp.resM = bsxfun(@minus,Dp.M,c);

        Dp1 = getrow(Dp,Dp.week==1);
        Dp4 = getrow(Dp,Dp.week==4);
        Dp12 = getrow(Dp,Dp.week==12);
        Dp24 = getrow(Dp,Dp.week==24);
        Dp52 = getrow(Dp,Dp.week==52);
        
        W1 = anova1(Dp1.resM);
        W4 = anova1(Dp4.resM);
        W12 = anova1(Dp12.resM);
        W24 = anova1(Dp24.resM);
        W52 = anova1(Dp52.resM);
        
%  % get individual patients and controls for each week
%         Dp  = tapply(Dp,{'subj_name','week'},{'M','mean(x,1)'});
%         Dc  = tapply(Dc,{'subj_name','week'},{'M','mean(x,1)'});
%         
%         getrow
%         
%         T = [];
%         
%         for p=unique(Dp.subj_name)'
%            S.W1 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==1&ismember(D.subj_name,p));
%            S.W4 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==4&ismember(D.subj_name,p));
%            S.W12 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==12&ismember(D.subj_name,p));
%            S.W24 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==24&ismember(D.subj_name,p));
%            S.W52 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==52&ismember(D.subj_name,p));         
%            
%         T = addstruct(T,S);
%         end
%         
%         for c=unique(Dc.subj_name)'          
%            V.C_W1 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==1&ismember(D.subj_name,c));
%            V.C_W4 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==4&ismember(D.subj_name,c));
%            V.C_W12 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==12&ismember(D.subj_name,c));
%            V.C_W24 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==24&ismember(D.subj_name,c));
%            V.C_W52 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==52&ismember(D.subj_name,c));         
%            
%         T = addstruct(T,V);
%         end
%         
% % take the residuals between patient/mean control & control/ mean control
% 
%      U = [];    
%       for i = 1:16
%           W1.PC_1 = corr(T.W1(i,:)',Dc_mean');
%           W1.delta_1 = T.W1(i,:)-Dc_mean;
%           U = addstruct(U,W1);
%       end
%       
%       for i = 1:17
%           W4.PC_4 = corr(T.W4(i,:)',Dc_mean');
%           W4.delta_4 = T.W4(i,:)-Dc_mean;
%           U = addstruct(U,W4);
%       end   
% 
%       for i = 1:17
%           W12.PC_12 = corr(T.W12(i,:)',Dc_mean');
%           W12.delta_12 = T.W12(i,:)-Dc_mean;
%           U = addstruct(U,W12);
%       end 
%       for i = 1:19
%           W24.PC_24 = corr(T.W24(i,:)',Dc_mean');
%           W24.delta_24 = T.W24(i,:)-Dc_mean;
%           U = addstruct(U,W24);
%       end
%       for i = 1:14
%           W52.PC_52 = corr(T.W52(i,:)',Dc_mean');
%           W52.delta_52 = T.W52(i,:)-Dc_mean;
%           U = addstruct(U,W52);
%       end
%        keyboard       
%         %build the averagefor each ROI pairing
%         
%       for i = 1:11
%           W.CC_1 = corr(T.C_W1(i,:)',Dc_mean');
%           W.C_delta_1 = T.C_W1(i,:)-Dc_mean;   
%           W.CC_4 = corr(T.C_W4(i,:)',Dc_mean');
%           W.C_delta_4 = T.C_W4(i,:)-Dc_mean;
%           W.CC_12 = corr(T.C_W12(i,:)',Dc_mean');
%           W.C_delta_12 = T.C_W12(i,:)-Dc_mean;
%           W.CC_24 = corr(T.C_W24(i,:)',Dc_mean');
%           W24.C_delta_24 = T.C_W24(i,:)-Dc_mean;
%           W.CC_52 = corr(T.C_W52(i,:)',Dc_mean');
%           W.C_delta_52 = T.C_W52(i,:)-Dc_mean;
%           U = addstruct(U,W);
%       end
%       
%       keyboard;
%       
%         
%       figure;
%         errorbar(mean(U.delta_1),std(U.delta_1))
%         hold on
%         errorbar(mean(U.C_delta_1),std(U.C_delta_1))
%         drawline(0,'dir','horz');
%       figure;
%         errorbar(mean(U.delta_4),std(U.delta_4))
%         hold on
%         errorbar(mean(U.C_delta_4),std(U.C_delta_4))
%         drawline(0,'dir','horz');
%       figure;
%         errorbar(mean(U.delta_12),std(U.delta_12))
%         hold on
%         errorbar(mean(U.C_delta_12),std(U.C_delta_12))
%         drawline(0,'dir','horz');
%       figure;
%         errorbar(mean(U.delta_24),std(U.delta_24))
%         hold on
%         errorbar(mean(U.C_delta_24),std(U.C_delta_24))
%         drawline(0,'dir','horz');
%       figure;
%         errorbar(mean(U.delta_52),std(U.delta_52))
%         hold on
%         errorbar(mean(U.C_delta_52),std(U.C_delta_52))
%         drawline(0,'dir','horz');
%         keyboard
%        
%         UU1 = mean(U.delta_1);
%         UU4 = mean(U.delta_4);
%         UU12 = mean(U.delta_12);
%         UU24 = mean(U.delta_24);
%         UU52 = mean(U.delta_52);
%         
%         T = [UU1', UU4', UU12',UU24',UU52'];
%         
%         % create Mixed-model for single ROIs?
     
    case 'split_mvc'
        %D   = load(fullfile(rootDir,'rs_networkCorr_lesionedIntra.mat'));                
        %D   = load(fullfile(rootDir,'rs_networkCorr_nonlesionedIntra.mat'));                
        D   = load(fullfile(ppDir,'rs_networkCorr.mat'));                
        
        
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
        
        % 0. Avg controls (all) and patients (per week and according to severity)
        C           = tapply(D,{'subj_name'},{'mat','mean(x,1)','subset',D.control==1});
                 
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
        keyboard;
    case 'split_mvc_hemis'
    fName = {'rs_networkCorr_lesionedIntra','rs_networkCorr_nonlesionedIntra'};
 
for i=1:length(fName)
        D = load(fullfile([ppDir],sprintf('%s.mat',fName{i})));
        
        for i=1:length(D.subj_name)
            m           = nonsym_squareform(D.M(i,:));
            D.mat(i,:)  = squareform(m);               %storing 10 pairs only
        end;
        
        %take only subgroup
        D   = getrow(D,(D.lesionType~=2)); %no cortical
        D   = getrow(D,(D.lesionType~=4)); %no mixed

        
        % 0. Avg controls (all) and patients (per week and according to severity)
        C           = tapply(D,{'subj_name'},{'mat','nanmean(x,1)','subset',D.control==1});
        [x,sn]      = pivottable(D.subj_name,D.week,D.FM,'nanmean(x)','subset',D.control==0);
        
        x           = nanmean(x(:,1:2),2);
        
        weeks = [1 4 12 24 52];
        
        for w = 1:length(weeks)
            
        Dp          = getrow(D,D.control==0 & D.week==weeks(w) & ismember(D.subj_name,sn(~isnan(x))));
        Dp.sev      = zeros(length(Dp.subj_name),1);
        
        Dp.sev(ismember(Dp.subj_name,sn(x<20)))            = 1;    % mvc <0.3 at weeks 1-4
        Dp.sev(ismember(Dp.subj_name,sn(x>=20)))           = 2;    % mvc >=0.3 at weeks 1-4
        
        P           = tapply(Dp,{'sev'},{'mat','nanmean(x,1)'},...
                               {'FM','nanmean(x)'});
                                    
        % 1. Compute correlations of
        %       - each control ens with left-out control ens mean
        %       - each control mm with left out control mm
        S = [];
        for i=1:length(C.subj_name)
            idxSubj     = i;
            idxLeftOut  = ~ismember(1:length(C.subj_name),i);

            rMat        = corr(C.mat(idxSubj,:)',nanmean(C.mat(idxLeftOut,:),1)');
            
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
        %ttest(S.rMat,T.rMat(T.sev==2),2,'paired');
        ttest(S.rMat,T.rMat(T.sev==2),2,'paired');
        fprintf('\n\n');
        end
        keyboard;
end        
    case 'split_mvc_hemis_ac'
         %before you run it, check: median or mean? MVC or FM
         %D      = load(fullfile(ppDir,'rs_networkCorr_nonlesionedIntra.mat'));
         D      = load(fullfile(ppDir,'rs_networkCorr_lesionedIntra.mat'));  
          
         %D       = load(fullfile(rootDir,'rs_networkCorr_lesionedIntra.mat'));                     
        
        for i=1:length(D.subj_name)
            m           = nonsym_squareform(D.M(i,:));
            D.mat(i,:)  = squareform(m);               %storing 10 pairs only
        end;
        
        %take only subgroup
        D   = getrow(D,(D.lesionType~=2)); %no cortical
        D   = getrow(D,(D.lesionType~=4)); %no mixed
        %D   = getrow(D,(D.lesionType~=3)); %no subcortical
        
        %   - mark weeks 1-12 as early, 24-52 as late
        week    = {[1,4],[24,52]};
        D.stage = zeros(length(D.subj_name),1);
        for i=1:2
            D.stage(ismember(D.week,week{i})) = i;
        end;
        D = getrow(D,D.stage ~= 0);
        

        % 0. Avg controls (all) and patients (per week and according to severity)
        C           = tapply(D,{'subj_name'},{'mat','nanmean(x,1)','subset',D.control==1});
        %C           = tapply(D,{'subj_name'},{'mat','nanmedian(x,1)','subset',D.control==1});
        
        %D.mvcNorm   = -D.mvcNorm;                       % for some reason mvcNorm is negative (to be fixed)
        %[x,sn]       = pivottable(D.subj_name,D.stage,D.mvcNorm,'nanmean(x)','subset',D.control==0);
        %[x,sn]      = pivottable(D.subj_name,D.stage,D.mvcNorm,'nanmedian(x)','subset',D.control==0);
        [x,sn]      = pivottable(D.subj_name,D.stage,D.FM,'nanmean(x)','subset',D.control==0);
        
        x            = x(:,1);
        
        stage = [1 2];
        
for w = 1:length(stage)
            
        Dp          = getrow(D,D.control==0 & D.stage==stage(w) & ismember(D.subj_name,sn(~isnan(x))));
        Dp.sev      = zeros(length(Dp.subj_name),1);
        
        %Dp.sev(ismember(Dp.subj_name,sn(x<0.3)))            = 1;    % mvc <0.3 at weeks 1-4
        %Dp.sev(ismember(Dp.subj_name,sn(x>=0.3)))           = 1;    % mvc >=0.3 at weeks 1-4
        
        Dp.sev(ismember(Dp.subj_name,sn(x<20)))            = 1;    % mvc <0.3 at weeks 1-4
        Dp.sev(ismember(Dp.subj_name,sn(x>=20)))           = 2;    % mvc >=0.3 at weeks 1-4
        
        %P           = tapply(Dp,{'sev'},{'mat','nanmean(x,1)'},...
        %           {'mvcNorm','nanmean(x)'});
        P           = tapply(Dp,{'sev'},{'mat','nanmean(x,1)'},...
                     {'FM','nanmean(x)'});                        
                     %{'mvcNorm','nanmedian(x)'});
                            
                                    
        % 1. Compute correlations of
        %       - each control ens with left-out control ens mean
        %       - each control mm with left out control mm
            S = [];
            for i=1:length(C.subj_name)
                idxSubj         = i;
                idxLeftOut      = ~ismember(1:length(C.subj_name),i);

                rMat            = corr(C.mat(idxSubj,:)',nanmean(C.mat(idxLeftOut,:),1)');

                Si.subj_name    = C.subj_name(i);
                Si.rMat         = rMat;
                Si.mat          = C.mat(idxSubj,:);
                S               = addstruct(S,Si);
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
            subplot (1,2,w)
            myboxplot(T.sev,T.rMat,'split',T.sev,'CAT',s(1).CAT,s(1).PLOT{:});
            drawline(mean(S.rMat),'dir','horz','color','b');
            drawline(mean(S.rMat)+stderr(S.rMat),'dir','horz','color','b','linestyle','--');
            drawline(mean(S.rMat)-stderr(S.rMat),'dir','horz','color','b','linestyle','--');
            title(stage(w))
            ylabel('pattern correlation')
            xlabel('severity')


            fprintf('Weeks %d\n',stage(w));
            fprintf('Sev 1\n');
            ttest(S.rMat,T.rMat(T.sev==1),2,'paired');
            fprintf('Sev 2\n');
            %ttest(S.rMat,T.rMat(T.sev==2),2,'paired');
            ttest(S.rMat,T.rMat(T.sev==2),2,'paired');
            fprintf('\n\n');
 end
            keyboard;
    case 'split_mvc_inter'
    %first check if median or mean and FM or mvc and split or no split
        
    D   = load(fullfile(ppDir,'rs_networkCorr.mat')); 
    
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
        C           = tapply(D,{'subj_name','week'},{'mat','nanmean(x,1)','subset',D.control==1});

        %D.mvcNorm   = -D.mvcNorm;                       % for some reason mvcNorm is negative (to be fixed)
        %[x,sn]      = pivottable(D.subj_name,D.week,D.mvcNorm,'nanmean(x)','subset',D.control==0);
        [x,sn]      = pivottable(D.subj_name,D.week,D.FM,'nanmean(x)','subset',D.control==0);
        
        x           = nanmean(x(:,1:2),2);
        
        weeks = [1 4 12 24 52];
        
        for w = 1:length(weeks)
            
        Dp          = getrow(D,D.control==0 & D.week==weeks(w) & ismember(D.subj_name,sn(~isnan(x))));
        Dp.sev      = zeros(length(Dp.subj_name),1);
        Dp.sev(ismember(Dp.subj_name,sn(x<20)))            = 1;    % mvc <0.3 at weeks 1-4
       %change accordingly
        Dp.sev(ismember(Dp.subj_name,sn(x>=20)))           = 2;    % mvc >=0.3 at weeks 1-4
        
        P           = tapply(Dp,{'sev','week'},{'mat','nanmedian(x,1)'},...
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
            Si.week         = C.week(idxSubj,:);
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
            Ti.week         = P.week;
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
    case 'split_mvc_inter_ac'
        %first check if median or mean and FM or mvc and split or no split      
    D   = load(fullfile(ppDir,'rs_networkCorr.mat'));  
             
    for i=1:length(D.subj_name)
            m           = nonsym_squareform(D.M(i,:));
            m           = m(:,1:5);
            m           = m(6:10,:);
            D.mat(i,:)  = nonsym_squareform(m);
    end;

        %take only subgroup
        D   = getrow(D,(D.lesionType~=2)); %no cortical
        D   = getrow(D,(D.lesionType~=4)); %no mixed
        
         %   - mark weeks 1-12 as early, 24-52 as late
        week    = {[1,4],[24,52]};
        D.stage = zeros(length(D.subj_name),1);
        for i=1:2
            D.stage(ismember(D.week,week{i})) = i;
        end;
        D = getrow(D,D.stage ~= 0);

        % 0. Avg controls (all) and patients (per week and according to severity)
        C            = tapply(D,{'subj_name'},{'mat','nanmean(x,1)','subset',D.control==1});

        %[x,sn]       = pivottable(D.subj_name,D.stage,D.mvcNorm,'nanmean(x)','subset',D.control==0);
        [x,sn]      = pivottable(D.subj_name,D.stage,D.FM,'nanmean(x)','subset',D.control==0);
        
        x           = x(:,1);
        
        stage = [1 2];
        
        for w = 1:length(stage)
            
        Dp          = getrow(D,D.control==0 & D.stage==stage(w) & ismember(D.subj_name,sn(~isnan(x))));
        Dp.sev      = zeros(length(Dp.subj_name),1);
        Dp.sev(ismember(Dp.subj_name,sn(x<20)))            = 1;    % mvc <0.3 at weeks 1-4
       %change accordingly
        Dp.sev(ismember(Dp.subj_name,sn(x>=20)))           = 2;    % mvc >=0.3 at weeks 1-4
        
        P           = tapply(Dp,{'sev'},{'mat','nanmean(x,1)'},...
                                        {'FM','nanmean(x)'});

        % 1. Compute correlations of
        %       - each control ens with left-out control ens mean
        %       - each control mm with left out control mm
        S = [];
        for i=1:length(C.subj_name)
            idxSubj     = i;
            idxLeftOut  = ~ismember(1:length(C.subj_name),i);

            rMat        = corr(C.mat(idxSubj,:)',nanmean(C.mat(idxLeftOut,:),1)');
            
            Si.subj_name    = C.subj_name(i);
            Si.rMat         = rMat;
            Si.mat          = C.mat(idxSubj,:);
            S               = addstruct(S,Si);
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
            T               = addstruct(T,Ti);
        end;
        varargout = {S,T};
        
        
        s   = style_sheet(sty_2grp_ac);
        subplot (1,2,w)
        myboxplot(T.sev,T.rMat,'split',T.sev,'CAT',s(1).CAT,s(1).PLOT{:}); 
        drawline(mean(S.rMat),'dir','horz','color','b');
        drawline(mean(S.rMat)+stderr(S.rMat),'dir','horz','color','b','linestyle','--');
        drawline(mean(S.rMat)-stderr(S.rMat),'dir','horz','color','b','linestyle','--');
        title(stage(w))
        ylabel('pattern correlation')
        %xlabel('severity') 
        
        fprintf('Weeks %d\n',stage(w));
        fprintf('Sev 1\n');
        ttest(S.rMat,T.rMat(T.sev==1),2,'paired');
         fprintf('Sev 2\n');
         ttest(S.rMat,T.rMat(T.sev==2),2,'paired');
         fprintf('\n\n');
        end 
    keyboard;   
    case 'split_mvcMagnitude'
        D      = load(fullfile(ppDir,'rs_networkCorr_nonlesionedIntra.mat'));             
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
        D       = load(fullfile(rootDir,'rs_networkCorr_nonlesionedIntra.mat'));        
        D       = rs_imana_110118('PP_removeInf',D);   
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
    
    case 'patients_changes_beh'
     D   = load(fullfile(ppDir,'nc_interhem.mat'));   
     D   = getrow(D,D.control==0);
         
     W1 = getrow(D,D.week==1);
     W4 = getrow(D,D.week==4);
     
     S = [];
 keyboard;    
    for i = unique(W1.subj_name)'
    [x y]= pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==4&ismember(D.subj_name,i));
    [Z]= pivottablerow(D.subj_name,D.FM,'nanmean','subset',D.week==4&ismember(D.subj_name,i));
    U.W4_M = x;
    U.subj_name = y;
    U.FM_4 = Z;
    
    S = addstruct(S,U);
    end 
     
    for j = unique(W4.subj_name)'
    [a b]= pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==1&ismember(D.subj_name,j));
    [ZZ]= pivottablerow(D.subj_name,D.FM,'nanmean','subset',D.week==1&ismember(D.subj_name,j));
    T.W1_M = a;
    T.subj_name_W1 = b;
    T.FM_1 = ZZ;
    S = addstruct(S,T);
    end 
    
    keyboard;
    
    Y = S.W1_M';
    Z = S.W4_M';
    X = mean(S.W4_M'-S.W1_M');
    
   
    YY = S.FM_1';
    ZZ = S.FM_4';
    XX = ZZ'-YY';
    scatterplot(XX,X','regression','linear','printcorr')
        
    case 'BEH_figure'
        D = load(fullfile(ppDir,'rs_preprocess.mat'));     
        
        % D = getrow(D,'subset of subcortical patients')
        T = tapply(D,{'subj_name','week','control'},{'FM','nanmean(x)'},...
                                                    {'ARAT','nanmean(x)'},...
                                                    {'ensOverall','nanmean(x)'},...
                                                    {'mvc','nanmean(x)'});
                                                
        % plot figures
        subplot(131);
        title('Recovery FM UE');
        lineplot(T.week,T.FM,'split',T.control,'plotfcn','nanmean');
        ylabel('FM score')
        xlabel('weeks')
        subplot(132);
        title('Recovery ARAT');
        lineplot(T.week,T.ARAT,'split',T.control,'plotfcn','nanmean');
        ylabel('Arat score')
        xlabel('weeks')
        subplot(133);
        title('Recovery hand strength');
        lineplot(T.week,T.mvc,'split',T.control,'plotfcn','nanmean');
        ylabel('hand strength in N')
        xlabel('weeks')
%         subplot(144);
%         lineplot(T.week,-T.ensOverall,'split',T.control,'plotfcn','nanmean');   
    case 'STATS_behavioural'
        D = rs_imana_110118('PP_studySubset');
        
        disp('FM');
        x=pivottable(D.subj_name,D.lesionType,D.FM,'nanmean','subset',D.week==1);
        ttest(x(:,1),x(:,2),2,'independent');
        
        disp('ARAT');
        x=pivottable(D.subj_name,D.lesionType,D.ARAT,'nanmean','subset',D.week==1);
        ttest(x(:,1),x(:,2),2,'independent');
        
        disp('MVC');
        x=pivottable(D.subj_name,D.lesionType,mean(D.mvc,2),'nanmean','subset',D.week==1);
        ttest(x(:,1),x(:,2),2,'independent');
    case 'STATS_behavioural_timecourse'
        D = rs_imana_110118('PP_studySubset');
        D.mvc = nanmean(D.mvc,2);
        
        D =     tapply(D,{'subj_name','week','control'},{'FM','nanmean'},...
                         {'ARAT','nanmean'},{'mvc','nanmean'});
        dsave(fullfile(statsDir,'behavioural_timecourse.dat'),D);       
    case 'STATS_timecourse_controls'
         D = load(fullfile(ppDir,'nc_interhem.mat'));        
         D = getrow(D,D.control==1);
         
         T = [];
         
        for  i = unique(D.subj_name)'
            W1 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==1&ismember(D.subj_name,i));
            W4 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==4&ismember(D.subj_name,i));
            W12 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==12&ismember(D.subj_name,i));
            W24 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==24&ismember(D.subj_name,i));
            W52 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==52&ismember(D.subj_name,i));
            
            
            if ~isempty(W1)&& ~isempty(W4)   
           Wa1_4 = corr(W1',W4');
           else
           Wa1_4 = NaN;
           end
           if ~isempty(W1)&& ~isempty(W12) 
           Wa1_12 = corr(W1',W12');
           else
           Wa1_12 = NaN;
           end
           if ~isempty(W1)&& ~isempty(W24) 
           Wa1_24 = corr(W1',W24');
           else
           Wa1_24 = NaN;
           end
           if ~isempty(W1)&& ~isempty(W52) 
           Wa1_52 = corr(W1',W52');
           else
           Wa1_52 = NaN;
           end
           
           if ~isempty(W1)&& ~isempty(W4) 
           W1_4 = corr(W1',W4');
           else
           W1_4 = NaN;
           end
           if~isempty(W4)&& ~isempty(W12) 
           W4_12 = corr(W4',W12');
           else
           W4_12 = NaN;
           end
           if ~isempty(W12)&& ~isempty(W24) 
           W12_24 = corr(W12',W24');
           else
           W12_24 = NaN;
           end
           if ~isempty(W24)&& ~isempty(W52) 
           W24_52 = corr(W24',W52');
           else
           W24_52 = NaN;
           end
           
           S.week_corr_acute    = [Wa1_4 Wa1_12 Wa1_24 Wa1_52];
           S.fz_week_acute      = fisherz(S.week_corr_acute);
           S.week_corr          = [W1_4 W4_12 W12_24 W24_52];
           S.fz_week_corr       = fisherz(S.week_corr);
           S.week               = [1 2 3 4];
           S.subj_name = i;
            
             T = addstruct(T,S);
        end  
        
 fprintf ('Inter') 
         save(fullfile(ppDir,'inter_timecourse.mat'),'-struct','T');   
         dsave(fullfile(statsDir,'inter_timecourse.dat'),T);  
% intra
        D = load(fullfile(ppDir,'nc_intrahem_lesioned.mat'));   %does not matter which one you load for controls     
        D = getrow(D,D.control==1);
        T = [];
        
       for  i = unique(D.subj_name)'
           W1 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==1&ismember(D.subj_name,i));
           W4 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==4&ismember(D.subj_name,i));
           W12 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==12&ismember(D.subj_name,i));
           W24 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==24&ismember(D.subj_name,i));
           W52 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==52&ismember(D.subj_name,i));
           
           if ~isempty(W1)&& ~isempty(W4)   
           Wa1_4 = corr(W1',W4');
           else
           Wa1_4 = NaN;
           end
           if ~isempty(W1)&& ~isempty(W12) 
           Wa1_12 = corr(W1',W12');
           else
           Wa1_12 = NaN;
           end
           if ~isempty(W1)&& ~isempty(W24) 
           Wa1_24 = corr(W1',W24');
           else
           Wa1_24 = NaN;
           end
           if ~isempty(W1)&& ~isempty(W52) 
           Wa1_52 = corr(W1',W52');
           else
           Wa1_52 = NaN;
           end
           
           if ~isempty(W1)&& ~isempty(W4) 
           W1_4 = corr(W1',W4');
           else
           W1_4 = NaN;
           end
           if~isempty(W4)&& ~isempty(W12) 
           W4_12 = corr(W4',W12');
           else
           W4_12 = NaN;
           end
           if ~isempty(W12)&& ~isempty(W24) 
           W12_24 = corr(W12',W24');
           else
           W12_24 = NaN;
           end
           if ~isempty(W24)&& ~isempty(W52) 
           W24_52 = corr(W24',W52');
           else
           W24_52 = NaN;
           end
           
           S.week_corr_acute    = [Wa1_4 Wa1_12 Wa1_24 Wa1_52];
           S.fz_week_acute      = fisherz(S.week_corr_acute);
           S.week_corr          = [W1_4 W4_12 W12_24 W24_52];
           S.fz_week_corr       = fisherz(S.week_corr);
           S.subj_name          = i;
           S.week               = [1 2 3 4];
           
            T = addstruct(T,S);
       end
       
 fprintf('Intra_lesioned');
  save(fullfile(ppDir,'intra_lesioned_timecourse.mat'),'-struct','T');  
  dsave(fullfile(statsDir,'intra_lesioned_timecourse.dat'),T); 
  
  % intra
        D = load(fullfile(ppDir,'nc_intrahem_nonlesioned.mat'));   %does not matter which one you load for controls     
        D = getrow(D,D.control==1);
        T = [];
        
       for  i = unique(D.subj_name)'
           W1 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==1&ismember(D.subj_name,i));
           W4 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==4&ismember(D.subj_name,i));
           W12 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==12&ismember(D.subj_name,i));
           W24 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==24&ismember(D.subj_name,i));
           W52 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==52&ismember(D.subj_name,i));
           
           
           if ~isempty(W1)&& ~isempty(W4)   
           Wa1_4 = corr(W1',W4');
           else
           Wa1_4 = NaN;
           end
           if ~isempty(W1)&& ~isempty(W12) 
           Wa1_12 = corr(W1',W12');
           else
           Wa1_12 = NaN;
           end
           if ~isempty(W1)&& ~isempty(W24) 
           Wa1_24 = corr(W1',W24');
           else
           Wa1_24 = NaN;
           end
           if ~isempty(W1)&& ~isempty(W52) 
           Wa1_52 = corr(W1',W52');
           else
           Wa1_52 = NaN;
           end
           
           if ~isempty(W1)&& ~isempty(W4) 
           W1_4 = corr(W1',W4');
           else
           W1_4 = NaN;
           end
           if~isempty(W4)&& ~isempty(W12) 
           W4_12 = corr(W4',W12');
           else
           W4_12 = NaN;
           end
           if ~isempty(W12)&& ~isempty(W24) 
           W12_24 = corr(W12',W24');
           else
           W12_24 = NaN;
           end
           if ~isempty(W24)&& ~isempty(W52) 
           W24_52 = corr(W24',W52');
           else
           W24_52 = NaN;
           end
           
           S.week_corr_acute    = [Wa1_4 Wa1_12 Wa1_24 Wa1_52];
           S.fz_week_acute      = fisherz(S.week_corr_acute);
           S.week_corr          = [W1_4 W4_12 W12_24 W24_52];
           S.fz_week_corr       = fisherz(S.week_corr);
           S.subj_name          = i;
           S.week               = [1 2 3 4];
           
            T = addstruct(T,S);
       end       
 fprintf('Intra_non_lesioned');
 save(fullfile(ppDir,'intra_nonlesioned_timecourse.mat'),'-struct','T'); 
 dsave(fullfile(statsDir,'intra_nonlesioned_timecourse.dat'),T); 
  
%all
        D = load(fullfile(ppDir,'nc_getThemall.mat'));   
        D = getrow(D,D.control==1);
        T = [];
        
       for  i = unique(D.subj_name)'
           W1 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==1&ismember(D.subj_name,i));
           W4 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==4&ismember(D.subj_name,i));
           W12 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==12&ismember(D.subj_name,i));
           W24 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==24&ismember(D.subj_name,i));
           W52 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==52&ismember(D.subj_name,i));
           
           if ~isempty(W1)&& ~isempty(W4)   
           Wa1_4 = corr(W1',W4');
           else
           Wa1_4 = NaN;
           end
           if ~isempty(W1)&& ~isempty(W12) 
           Wa1_12 = corr(W1',W12');
           else
           Wa1_12 = NaN;
           end
           if ~isempty(W1)&& ~isempty(W24) 
           Wa1_24 = corr(W1',W24');
           else
           Wa1_24 = NaN;
           end
           if ~isempty(W1)&& ~isempty(W52) 
           Wa1_52 = corr(W1',W52');
           else
           Wa1_52 = NaN;
           end
           
           if ~isempty(W1)&& ~isempty(W4) 
           W1_4 = corr(W1',W4');
           else
           W1_4 = NaN;
           end
           if~isempty(W4)&& ~isempty(W12) 
           W4_12 = corr(W4',W12');
           else
           W4_12 = NaN;
           end
           if ~isempty(W12)&& ~isempty(W24) 
           W12_24 = corr(W12',W24');
           else
           W12_24 = NaN;
           end
           if ~isempty(W24)&& ~isempty(W52) 
           W24_52 = corr(W24',W52');
           else
           W24_52 = NaN;
           end
           
           S.week_corr_acute    = [Wa1_4 Wa1_12 Wa1_24 Wa1_52];
           S.fz_week_acute      = fisherz(S.week_corr_acute);
           S.week_corr          = [W1_4 W4_12 W12_24 W24_52];
           S.fz_week_corr       = fisherz(S.week_corr);
           S.subj_name          = i;
           S.week               = [1 2 3 4];
           
            T = addstruct(T,S);
       end       
 fprintf('all');  
 save(fullfile(ppDir,'all_timecourse.mat'),'-struct','T'); 
 dsave(fullfile(statsDir,'all_timecourse.dat'),T);
 
 keyboard;
    case 'STATS_timecourse_patients'    
         D = load(fullfile(ppDir,'nc_interhem.mat'));        
         D = getrow(D,D.control==0);
         T = [];
         
        for  i = unique(D.subj_name)'
            W1 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==1&ismember(D.subj_name,i));
            W4 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==4&ismember(D.subj_name,i));
            W12 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==12&ismember(D.subj_name,i));
            W24 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==24&ismember(D.subj_name,i));
            W52 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==52&ismember(D.subj_name,i));
            
            
            if ~isempty(W1)&& ~isempty(W4)   
           Wa1_4 = corr(W1',W4');
           else
           Wa1_4 = NaN;
           end
           if ~isempty(W1)&& ~isempty(W12) 
           Wa1_12 = corr(W1',W12');
           else
           Wa1_12 = NaN;
           end
           if ~isempty(W1)&& ~isempty(W24) 
           Wa1_24 = corr(W1',W24');
           else
           Wa1_24 = NaN;
           end
           if ~isempty(W1)&& ~isempty(W52) 
           Wa1_52 = corr(W1',W52');
           else
           Wa1_52 = NaN;
           end
           
           if ~isempty(W1)&& ~isempty(W4) 
           W1_4 = corr(W1',W4');
           else
           W1_4 = NaN;
           end
           if~isempty(W4)&& ~isempty(W12) 
           W4_12 = corr(W4',W12');
           else
           W4_12 = NaN;
           end
           if ~isempty(W12)&& ~isempty(W24) 
           W12_24 = corr(W12',W24');
           else
           W12_24 = NaN;
           end
           if ~isempty(W24)&& ~isempty(W52) 
           W24_52 = corr(W24',W52');
           else
           W24_52 = NaN;
           end
           
           S.week_corr_acute    = [Wa1_4 Wa1_12 Wa1_24 Wa1_52];
           S.fz_week_acute      = fisherz(S.week_corr_acute);
           S.week_corr          = [W1_4 W4_12 W12_24 W24_52];
           S.fz_week_corr       = fisherz(S.week_corr);
           S.subj_name          = i;
           S.week               = [1 2 3 4];
            
             T = addstruct(T,S);
        end  
        
 fprintf ('Inter') 
         save(fullfile(ppDir,'inter_timecourse_patient.mat'),'-struct','T');   
         dsave(fullfile(statsDir,'inter_timecourse_patient.dat'),T);  
% intra
        D = load(fullfile(ppDir,'nc_intrahem_lesioned.mat'));   %does not matter which one you load for controls     
        D = getrow(D,D.control==0);
        T = [];
        
       for  i = unique(D.subj_name)'
           W1 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==1&ismember(D.subj_name,i));
           W4 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==4&ismember(D.subj_name,i));
           W12 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==12&ismember(D.subj_name,i));
           W24 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==24&ismember(D.subj_name,i));
           W52 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==52&ismember(D.subj_name,i));
           
           if ~isempty(W1)&& ~isempty(W4)   
           Wa1_4 = corr(W1',W4');
           else
           Wa1_4 = NaN;
           end
           if ~isempty(W1)&& ~isempty(W12) 
           Wa1_12 = corr(W1',W12');
           else
           Wa1_12 = NaN;
           end
           if ~isempty(W1)&& ~isempty(W24) 
           Wa1_24 = corr(W1',W24');
           else
           Wa1_24 = NaN;
           end
           if ~isempty(W1)&& ~isempty(W52) 
           Wa1_52 = corr(W1',W52');
           else
           Wa1_52 = NaN;
           end
           
           if ~isempty(W1)&& ~isempty(W4) 
           W1_4 = corr(W1',W4');
           else
           W1_4 = NaN;
           end
           if~isempty(W4)&& ~isempty(W12) 
           W4_12 = corr(W4',W12');
           else
           W4_12 = NaN;
           end
           if ~isempty(W12)&& ~isempty(W24) 
           W12_24 = corr(W12',W24');
           else
           W12_24 = NaN;
           end
           if ~isempty(W24)&& ~isempty(W52) 
           W24_52 = corr(W24',W52');
           else
           W24_52 = NaN;
           end
           
           S.week_corr_acute    = [Wa1_4 Wa1_12 Wa1_24 Wa1_52];
           S.fz_week_acute      = fisherz(S.week_corr_acute);
           S.week_corr          = [W1_4 W4_12 W12_24 W24_52];
           S.fz_week_corr       = fisherz(S.week_corr);
           S.subj_name          = i;
           S.week               = [1 2 3 4];
           
            T = addstruct(T,S);
       end
       
 fprintf('Intra_lesioned');
  save(fullfile(ppDir,'intra_lesioned_timecourse_patient.mat'),'-struct','T');  
  dsave(fullfile(statsDir,'intra_lesioned_timecourse_patient.dat'),T); 
  
  % intra
        D = load(fullfile(ppDir,'nc_intrahem_nonlesioned.mat'));   %does not matter which one you load for controls     
        D = getrow(D,D.control==0);
        T = [];
        
       for  i = unique(D.subj_name)'
           W1 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==1&ismember(D.subj_name,i));
           W4 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==4&ismember(D.subj_name,i));
           W12 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==12&ismember(D.subj_name,i));
           W24 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==24&ismember(D.subj_name,i));
           W52 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==52&ismember(D.subj_name,i));
           
           
           if ~isempty(W1)&& ~isempty(W4)   
           Wa1_4 = corr(W1',W4');
           else
           Wa1_4 = NaN;
           end
           if ~isempty(W1)&& ~isempty(W12) 
           Wa1_12 = corr(W1',W12');
           else
           Wa1_12 = NaN;
           end
           if ~isempty(W1)&& ~isempty(W24) 
           Wa1_24 = corr(W1',W24');
           else
           Wa1_24 = NaN;
           end
           if ~isempty(W1)&& ~isempty(W52) 
           Wa1_52 = corr(W1',W52');
           else
           Wa1_52 = NaN;
           end
           
           if ~isempty(W1)&& ~isempty(W4) 
           W1_4 = corr(W1',W4');
           else
           W1_4 = NaN;
           end
           if~isempty(W4)&& ~isempty(W12) 
           W4_12 = corr(W4',W12');
           else
           W4_12 = NaN;
           end
           if ~isempty(W12)&& ~isempty(W24) 
           W12_24 = corr(W12',W24');
           else
           W12_24 = NaN;
           end
           if ~isempty(W24)&& ~isempty(W52) 
           W24_52 = corr(W24',W52');
           else
           W24_52 = NaN;
           end
           
           S.week_corr_acute    = [Wa1_4 Wa1_12 Wa1_24 Wa1_52];
           S.fz_week_acute      = fisherz(S.week_corr_acute);
           S.week_corr          = [W1_4 W4_12 W12_24 W24_52];
           S.fz_week_corr       = fisherz(S.week_corr);
           S.subj_name          = i;
           S.week               = [1 2 3 4];
           
            T = addstruct(T,S);
       end       
 fprintf('Intra_non_lesioned');
 save(fullfile(ppDir,'intra_nonlesioned_timecourse_patient.mat'),'-struct','T'); 
 dsave(fullfile(statsDir,'intra_nonlesioned_timecourse_patient.dat'),T); 
  
%all
        D = load(fullfile(ppDir,'nc_getThemall.mat'));   
        D = getrow(D,D.control==0);
        T = [];
        
       for  i = unique(D.subj_name)'
           W1 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==1&ismember(D.subj_name,i));
           W4 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==4&ismember(D.subj_name,i));
           W12 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==12&ismember(D.subj_name,i));
           W24 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==24&ismember(D.subj_name,i));
           W52 = pivottablerow(D.subj_name,D.M,'nanmean(x,1)','subset',D.week==52&ismember(D.subj_name,i));
           
           if ~isempty(W1)&& ~isempty(W4)   
           Wa1_4 = corr(W1',W4');
           else
           Wa1_4 = NaN;
           end
           if ~isempty(W1)&& ~isempty(W12) 
           Wa1_12 = corr(W1',W12');
           else
           Wa1_12 = NaN;
           end
           if ~isempty(W1)&& ~isempty(W24) 
           Wa1_24 = corr(W1',W24');
           else
           Wa1_24 = NaN;
           end
           if ~isempty(W1)&& ~isempty(W52) 
           Wa1_52 = corr(W1',W52');
           else
           Wa1_52 = NaN;
           end
           
           if ~isempty(W1)&& ~isempty(W4) 
           W1_4 = corr(W1',W4');
           else
           W1_4 = NaN;
           end
           if~isempty(W4)&& ~isempty(W12) 
           W4_12 = corr(W4',W12');
           else
           W4_12 = NaN;
           end
           if ~isempty(W12)&& ~isempty(W24) 
           W12_24 = corr(W12',W24');
           else
           W12_24 = NaN;
           end
           if ~isempty(W24)&& ~isempty(W52) 
           W24_52 = corr(W24',W52');
           else
           W24_52 = NaN;
           end
           
           S.week_corr_acute    = [Wa1_4 Wa1_12 Wa1_24 Wa1_52];
           S.fz_week_acute      = fisherz(S.week_corr_acute);
           S.week_corr          = [W1_4 W4_12 W12_24 W24_52];
           S.fz_week_corr       = fisherz(S.week_corr);
           S.subj_name          = i;
           S.week               = [1 2 3 4];
           
            T = addstruct(T,S);
       end       
 fprintf('all');  
 save(fullfile(ppDir,'all_timecourse_patient.mat'),'-struct','T'); 
 dsave(fullfile(statsDir,'all_timecourse_patient.dat'),T);
 
 keyboard;
    case 'timecourse_FIG'
    A = load(fullfile(ppDir,'all_timecourse.mat'));   
    AA = (A.week_corr_acute);
    r = nanmean(fisherz(AA),2);
    m = nanmean(r); SE = stderr(r); 
    fprintf('ALL_timecourse_controls r = %2.3f (%2.3f - %2.3f)\n',(m),(m-1.96*SE),(m+1.96*SE));
    
    B = load(fullfile(ppDir,'all_timecourse_patient.mat'));   
    BB = (B.week_corr_acute);
    r = nanmean(fisherz(BB),2);
    m = nanmean(r); SE = stderr(r); 
    fprintf('ALL_timecourse_patients r = %2.3f (%2.3f - %2.3f)\n',(m),(m-1.96*SE),(m+1.96*SE));
    
    C = load(fullfile(ppDir,'inter_timecourse.mat'));   
    CC = (C.week_corr_acute);
    r = nanmean(fisherz(CC),2);
    m = nanmean(r); SE = stderr(r); 
    fprintf('INTER_timecourse_controls r = %2.3f (%2.3f - %2.3f)\n',(m),(m-1.96*SE),(m+1.96*SE));
    
    D = load(fullfile(ppDir,'inter_timecourse_patient.mat'));   
    DD = (D.week_corr_acute);
    r = nanmean(fisherz(DD),2);
    m = nanmean(r); SE = stderr(r); 
    fprintf('INTER_timecourse_patients r = %2.3f (%2.3f - %2.3f)\n',(m),(m-1.96*SE),(m+1.96*SE));
    
    E = load(fullfile(ppDir,'intra_lesioned_timecourse.mat'));   
    EE = (E.week_corr_acute);
    r = nanmean(fisherz(EE),2);
    m = nanmean(r); SE = stderr(r); 
    fprintf('INTRA_timecourse_controls r = %2.3f (%2.3f - %2.3f)\n',(m),(m-1.96*SE),(m+1.96*SE));
    
    F = load(fullfile(ppDir,'intra_lesioned_timecourse_patient.mat'));   
    FF = (F.week_corr_acute);
    r = nanmean(fisherz(FF),2);
    m = nanmean(r); SE = stderr(r); 
    fprintf('INTRA_timecourse_patients r = %2.3f (%2.3f - %2.3f)\n',(m),(m-1.96*SE),(m+1.96*SE));
    
    G = load(fullfile(ppDir,'intra_nonlesioned_timecourse.mat'));   
    GG = (G.week_corr_acute);
    r = nanmean(fisherz(GG),2);
    m = nanmean(r); SE = stderr(r); 
    fprintf('INTRA_NON_timecourse_controls r = %2.3f (%2.3f - %2.3f)\n',(m),(m-1.96*SE),(m+1.96*SE));
    
    H = load(fullfile(ppDir,'intra_nonlesioned_timecourse_patient.mat'));   
    HH = (H.week_corr_acute);
    r = nanmean(fisherz(HH),2);
    m = nanmean(r); SE = stderr(r); 
    fprintf('INTRA_NON_timecourse_patients r = %2.3f (%2.3f - %2.3f)\n',(m),(m-1.96*SE),(m+1.96*SE));
    
    
    keyboard;
    
    case 'STATS_controlInterSubj'
        D = load(fullfile(ppDir,'nc_interhem.mat')); 
        D = getrow(D,D.control==1);
        x = pivottablerow(D.subj_name,D.M,'nanmean(x,1)');
        r = corr(x') + diag(nan(1,11));
        r = nanmean(fisherz(r),2);
        m = mean(r); SE = stderr(r); 
        fprintf('Interhem r = %2.3f (%2.3f - %2.3f)\n',fisherinv(m),fisherinv(m-1.96*SE),fisherinv(m+1.96*SE));
        
        D = load(fullfile(ppDir,'nc_intrahem_lesioned.mat')); 
        D = getrow(D,D.control==1);
        x = pivottablerow(D.subj_name,D.M,'nanmean(x,1)');
        r = corr(x') + diag(nan(1,11));
        r = nanmean(fisherz(r),2);
        m = mean(r); SE = stderr(r);
        fprintf('Intra-lesioned r = %2.3f (%2.3f - %2.3f)\n',fisherinv(m),fisherinv(m-1.96*SE),fisherinv(m+1.96*SE));
    case 'STATS_split_half'
        
         S = load(fullfile(ppDir,'nc_interhem_splithalf.mat'));        
         S = getrow(S,S.lesionType ~=2);
         S = getrow(S,S.lesionType ~=4);
         [x,y]=pivottable([S.subj_name],[],S.sh_fz_r,'(nanmean(x))','subset',S.control==0);
         [z,y]=pivottable([S.subj_name],[],S.sh_fz_r,'(nanmean(x))','subset',S.control==1);
         

        mp = nanmean(x); SEp = stderr(x);
        mc = nanmean(z); SEc = stderr(z); 
        fprintf('Interhem patient r = %2.3f (%2.3f - %2.3f)\n',fisherinv(mp),fisherinv(mp-1.96*SEp),fisherinv(mp+1.96*SEp));
        fprintf('Interhem control r = %2.3f (%2.3f - %2.3f)\n',fisherinv(mc),fisherinv(mc-1.96*SEc),fisherinv(mc+1.96*SEc));
        
        S = load(fullfile(ppDir,'nc_intrahem_lesioned_splithalf.mat')); 
        S = getrow(S,S.lesionType ~=2);
        S = getrow(S,S.lesionType ~=4);
        [x,y]=pivottable([S.subj_name],[],S.sh_fz_r,'(nanmean(x))','subset',S.control==0);
        [z,y]=pivottable([S.subj_name],[],S.sh_fz_r,'(nanmean(x))','subset',S.control==1);
        
        
        mp = nanmean(x); SEp = stderr(x);
        mc = nanmean(z); SEc = stderr(z); 
        fprintf('Intrahem_les patient r = %2.3f (%2.3f - %2.3f)\n',fisherinv(mp),fisherinv(mp-1.96*SEp),fisherinv(mp+1.96*SEp));
        fprintf('Intrahem_les control r = %2.3f (%2.3f - %2.3f)\n',fisherinv(mc),fisherinv(mc-1.96*SEc),fisherinv(mc+1.96*SEc));
       
        S = load(fullfile(ppDir,'nc_intrahem_nonlesioned_splithalf.mat')); 
        S = getrow(S,S.lesionType ~=2);
        S = getrow(S,S.lesionType ~=4);
        [x,y]=pivottable([S.subj_name],[],S.sh_fz_r,'(nanmean(x))','subset',S.control==0);
        [z,y]=pivottable([S.subj_name],[],S.sh_fz_r,'(nanmean(x))','subset',S.control==1);
        
        mp = nanmean(x); SEp = stderr(x);
        mc = nanmean(z); SEc = stderr(z); 
        fprintf('Intrahem_nonles patient r = %2.3f (%2.3f - %2.3f)\n',fisherinv(mp),fisherinv(mp-1.96*SEp),fisherinv(mp+1.96*SEp));
        fprintf('Intrahem_nonles control r = %2.3f (%2.3f - %2.3f)\n',fisherinv(mc),fisherinv(mc-1.96*SEc),fisherinv(mc+1.96*SEc));     
        
    case 'STATS_make_dunnetts_testfile'
        C   = dload('MainFigure_all_EUC_C.dat');
        P   = dload('MainFigure_all_EUC_P.dat');
        D   = addstruct(C,P);
        dsave('MainFigure_all_EUC.dat',D);
        keyboard;
        
    case 'otherwise'
        disp('No such case...');
     
    %case 'PP_getSubset'                 % get subset of data from the main data structure
%         exclROI     = [6:8];                % regions to exclude (no parietal or V1)
%         exclSubj    = {''};                 % remove these participants
%         type        = 'lesioned';           % lesioned/nonlesioned/inter/all
%         fname       = 'rs_preprocess.mat';  % over-ride file name to open a different structure
%         lesionType  = [];                   % all lesion types
%         vararginoptions(varargin,{'exclROI','exclSubj','type','fname','lesionType'});
%         
%         % 1. Load data for split-half reliabilities
%         %       - only include required ROIs
%         D   = load(fullfile(ppDir,fname));
%         idx = ~(ismember(D.regType(:,1),exclROI) | ismember(D.regType(:,2),exclROI));
%         D   = getrow(D,idx);
%         
%         %       - remove these subjects
%         idx = ~ismember(D.subj_name,exclSubj);
%         D   = getrow(D,idx);
% 
%         %       - get rows corresponding only to the type
%         switch(type)
%             case 'lesioned'
%                 Dc  = getrow(D,D.control==1 & D.hl~=0);                 % controls, intra-cortical only
%                 Dp  = getrow(D,D.control==0 & D.lesionSide==D.hl);      % patients, lesioned hem, intra-cortical only.         
%             case 'nonlesioned'
%                 Dc  = getrow(D,D.control==1 & D.hl~=0);                 % controls, intra-cortical only
%                 Dp  = getrow(D,D.control==0 & D.lesionSide~=D.hl);      % patients, nonlesioned hem, intra-cortical only.         
%             case 'inter'
%                 Dc  = getrow(D,D.control==1 & D.hl==0);                 % controls, inter-cortical only
%                 Dp  = getrow(D,D.control==0 & D.hl==0);                 % patients, inter-cortical only                 
%             case 'all'
%                 Dc  = getrow(D,D.control==1);                           % all controls
%                 Dp  = getrow(D,D.control==0);                           % all patients
%         end;
%         
%         %       - only get specified lesion types
%         if ~isempty(lesionType)
%             Dp  = getrow(Dp,ismember(Dp.lesionType,lesionType));
%         end;
%         
%         % 2. return structure
%         T   = addstruct(Dc,Dp);        
%         varargout = {T};
        
    %case 'PP_removeSubj'                % some participants have bad correlation estimates for multiple reasons
%         D = varargin{1};
%         
%         % 0. Figure out which participants have inf in their fisherz scores
%         %   - this usually happens when two rois have no data,
%         %   region_getdata returns 0's or 1's for all voxel values, leading
%         %   to a correlation r=1, and a corresponding fisherz r=Inf,-Inf
%         [x,sn,week]     = pivottable(D.subj_name,D.week,D.fz_r,'sum(isinf(x))');
%         [i,j]           = ind2sub(size(x),find(x>0));
%         
%         idx = zeros(length(D.subj_name),1);    
%         for k=1:length(i)
%             l       = strcmp(D.subj_name,sn(i(k))) & D.week==week(j(k));
%             idx(l)  = 1;
%             fprintf('%s, week %d removed\n',sn{i(k)},week(j(k)));
%         end;
%         D = getrow(D,~idx);
%         varargout = {D};    
    %case 'PP_additional'        % additional variables added to preprocessed data to make analysis easier
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