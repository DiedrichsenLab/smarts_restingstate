function varargout=rs_imana(what,varargin)
% MB - analysis script for the SMARTS resting state data


rootDir               = '/Volumes/porsche/data/smarts';
rootDir               = '/Volumes/MotorControl/data/smarts';
rootDir               = '/Volumes/diedrichsen_data$/data/smarts';
behDir              = [rootDir '/bedside/analysis'];
rsDir               = [rootDir '/fmri/restingstate_imaging'];
compareDir          = [rootDir '/fmri/restingstate_imaging/preprocess/files'];
roiDir              = [rootDir '/fmri/RegionOfInterest'];
ppDir               = [rsDir '/preprocess'];
analysisDir         = [rsDir '/preprocess'];
preZipDir           = [rootDir '/fmri/restingstate_imaging/pre'];
statsDir            = [rootDir '/fmri/restingstate_imaging/preprocess/R/'];
roiPark2011         = [rsDir '/RegionOfInterest-Park2011'];

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

lesionType = {'none','cortical','subcortical','mixed'}; %1,2,3,4
HemN = {'LeftHem','RightHem'};
hem = {'lh','rh'}; 
        
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
    case 'preprocess_comparison'              % Compares preprocessing
        % 0. check which subjects have the raw and preprocessed files
        D = rs_preprocess_ana('check_data_exists');
        D = getrow(D,D.existPPImg & D.existRawImg);
        D = getrow(D,D.lesionside==0);                                      % get only controls
        
        idx = ismember(D.subj_name,'CUP_1002') & D.week==24;
        D   = getrow(D,~idx);
        
        % need to estimate rel. over the two preprocessing techniques
        %   pre:        ICA-based noise removal pipeline
        %   preproc:    regression pipeline
        type = {'pre_Bold_Rest','preproc_Bold_Rest'};
        %roi  = [2 10];  % left and right M1                                 % change if you only want M1
        
        S = [];
        for i=1:length(D.SN)
            Di = getrow(D,i);                                               % get i-th session
            fprintf('%s - %s\n',Di.subj_name{1},Di.Week{1});
            
            % load ROI definitions for that subjection
            R   = load(fullfile(roiDir,sprintf('%s_regions_all.mat',Di.subj_name{1})));
            
            % loop over two preprocessing types
            for p=1:length(type)
                %fdir    = fullfile(rsDir,Di.subj_name{1},Di.Week{1});
                fdir    = fullfile(compareDir);
                fname   = sprintf('%s_%s_%s.nii',type{p},Di.subj_name{1},Di.Week{1});
                
                if exist(fullfile(fdir,fname),'file')
                    v   = spm_vol(fullfile(fdir,fname));                    % load volumes from rs data
                    %ts = region_getdata(v,R.R(roi));                       % get time series for M1-left and M1-right
                    ts = region_getdata(v,R.R);
                    % calculate corr for time series
                    ts1 = mean(ts{1},2);
                    ts2 = mean(ts{2},2);
                    r   = corr(ts1,ts2);
                else
                    fprintf('file missing\n');
                    r   = nan;
                end;
                
                % write out data to structure
                Si              = Di;
                Si.ptype        = p;
                Si.pname        = type(p);
                Si.r            = r;
                S               = addstruct(S,Si);
            end;
        end;
        varargout = {S};
        %save(fullfile(analysisDir,'rs_compare_preprocess.mat'),'-struct','S');
        save(fullfile(analysisDir,'rs_compare_preprocess_whole_pattern.mat'),'-struct','S');
    case 'FIG_compare_prepro'
        D = load(fullfile(ppDir,'rs_compare_preprocess_whole_pattern.mat'));
        %
        s   = style_sheet(sty_2grp_ac,'leg',{'preproc_Bold_Rest', 'pre_Bold_Rest'},'leglocation','northoutside');
        h = make_figure;
        %   subplot(3,2,1)
        title('Comparison noise correction')
        lineplot(D.week,D.r,'plotfcn','nanmean','errorfcn','stderr','split',D.ptype,'CAT',s(1).CAT,s(1).PLOT{:});ylabel('correlation index')
        xlabel('weeks')
        dsave(fullfile(statsDir,'rs_compare_preprocess_whole_pattern_stat.dat'),D);
        keyboard;
        
        x = pivottable([D.subj_name],[D.ptype],D.r,'nanmean');
        preproc      = nanmean(x(:,1)); SEp = stderr(x(:,1));
        pre          = nanmean(x(:,2)); SEc = stderr(x(:,2));
        
        fprintf('Preproc r = %2.3f (%2.3f - %2.3f)\n',fisherinv(preproc),fisherinv(preproc-1.96*SEp),fisherinv(preproc+1.96*SEp));
        fprintf('Pre r = %2.3f (%2.3f - %2.3f)\n',fisherinv(pre),fisherinv(pre-1.96*SEc),fisherinv(pre+1.96*SEc));
        
    case 'PP_rawTimeSeries'                   % Extraction of fisher-z correlations from different regions of interest from "freesurfer atlas"
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
                %fDir= fullfile(rsDir,D.subj_name{i},D.Week{i});
                %v   = spm_vol(fullfile(fDir,fname));       % load volumes from rs data
                v   = spm_vol(fullfile(compareDir,fname));       % load volumes from rs data
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
                    Si.regType  = regType(Si.rpair)';
                    
                    S           = addstruct(S,Si);
                end;
            end;
        end;
        varargout = {S};
        save(fullfile(ppDir,'rs_preprocess_test.mat'),'-struct','S');
    case 'PP_rawTimeSeries_splitHalf'
        
        pre = {'preproc_Bold_Rest'};
        
        lesionType = {'none','cortical','subcortical','mixed'}; %1,2,3,4
        roi  = combnk(1:16,2); % all pairs
        
        % 0. check which subjects have the raw and preprocessed files
        D = rs_preprocess_ana('check_data_exists');
        D = getrow(D,D.existPPImg & D.existRawImg);
        D = getrow(D,D.lesionside==0);
        
        % remove bad control time point
        idx     = ismember(D.subj_name,'CUP_1002') & D.week==24;
        D      = getrow(D,~idx);
        
        %       D   = dload(fullfile(ppDir,'patient_list.txt'));
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
                %fDir= fullfile(rsDir,D.subj_name{i},D.Week{i});
                fDir = compareDir;
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
                Si.lesionSide   = D.lesionside(i);                          %lesionside = 1 = left
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
    case 'PP_voxelwise_define_region' 
        % 0. check which subjects have the raw and preprocessed files
        D = rs_preprocess_ana('check_data_exists');
        D = getrow(D,D.existPPImg & D.existRawImg);
        [names,ia,ib] = unique(D.subj_name); 
        
        % Define the area of the damaged hemisphere that we will consider
        atlasSurfDir = [rootDir '/fmri/surface_caret/fsaverage_sym/LeftHem'];
        FL = caret_load(fullfile(atlasSurfDir,'lh.FLAT.coord'));
        nodes = find(FL.data(:,1)>-25 & FL.data(:,1)<18 & FL.data(:,2)<42 & FL.data(:,2)>-27);
        
        for i = ia' 
            fprintf('defining region for %s \n',D.subj_name{i}); 
            regFile = fullfile(roiDir,sprintf('%s_regions_all.mat',D.subj_name{i}));
            R   = load(fullfile(roiDir,sprintf('%s_regions_all.mat',D.subj_name{i}))); % load subject region definition
            s = regexp(R.R{1}.image, '/', 'split');
            mask = fullfile(rootDir,'fmri','anatomicals',D.subj_name{i},s{end-1},s{end});
                
            % Define the healthy and lesioned hemisphere
            sideL = D.lesionside(i);
            if (sideL == 0)
                sideL = 1;
            end;
            sideN = 3-sideL;
                
            % Make the region structure
            N.type = 'surf_nodes'; 
            N.location = nodes; 
            N.white = fullfile(rootDir,'fmri','surface_caret',['x' D.subj_name{i}],HemN{sideL},[hem{sideL} '.WHITE.coord']);
            N.pial = fullfile(rootDir,'fmri','surface_caret',['x' D.subj_name{i}],HemN{sideL},[hem{sideL} '.PIAL.coord']);
            N.image = mask; 
            N.name = [D.subj_name{i} 'voxelwise']; 
            N.linedef = [5 0 1];
            
            % Caluculate and save 
            N = region_calcregions(N); 
            name = fullfile(roiDir,sprintf('%s_regions_voxelwise.mat',D.subj_name{i})); % load subject region definition
            save(name,'-struct','N'); 
        end; 
            
    case 'PP_voxelWiseCorrelation'                      % Extraction of fisher-z correlations from different regions of interest from "freesurfer atlas"
        pre = {'preproc_Bold_Rest'};

        % 0. check which subjects have the raw and preprocessed files
        D = rs_preprocess_ana('check_data_exists');
        D = getrow(D,D.existPPImg & D.existRawImg);
        
        %         D   = dload(fullfile(rsDir,'patient_list.txt'));
        B   = load(fullfile(behDir,'ens_alldat_mvc.mat'));
        SS  = dload(fullfile(rootDir,'bedside','data_clean','subject_list.txt'));
        L   = dload(fullfile(rootDir,'DTI','lesion_analysis','lesion.dat'));
        
        SS.subj_name    = strcat(SS.Centre,'_',num2str(SS.ID));
        L.subj_name     = strcat(L.Centre,'_',num2str(L.ID));
        
        % Define the area of the damaged hemisphere that we will consider
        atlasSurfDir = [rootDir '/fmri/surface_caret/fsaverage_sym/LeftHem'];
        FL = caret_load(fullfile(atlasSurfDir,'lh.FLAT.coord'));
        nodes = FL.data(:,1)>-25 & FL.data(:,1)<18 & FL.data(:,2)<42 & FL.data(:,2)>-27;
        
        S = [];
        for i=1:length(D.SN)
            fname   = strcat(pre{1},'_',D.Centre{i},'_',num2str(D.ID(i)),'_',D.Week{i},'.nii');
            sn      = strcat(D.Centre{i},'_',num2str(D.ID(i)));
            fprintf('%s\n',fname);
            
            regFile = fullfile(roiDir,sprintf('%s_regions_all.mat',sn));
            if ~exist(regFile)
                fprintf('No region definition for this subject...\n');
            else
                R   = load(fullfile(roiDir,sprintf('%s_regions_all.mat',sn))); % load subject region definition
                N   = load(fullfile(roiDir,sprintf('%s_regions_voxelwise.mat',sn))); % load subject region definition
                
                % Define the healthy and lesioned hemisphere
                sideL = D.lesionside(i);
                if (sideL == 0)
                    sideL = 1;
                end;
                sideN = 3-sideL;
                regions = {R.R{2+(sideN-1)*8},N};
                % Extract data 
                v   = spm_vol(fullfile(compareDir,fname));       % load volumes from rs data
                ts  = region_getdata(v,regions);                % get time series for all ROIs
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
                Si.lesionSide   = D.lesionside(i);                          %lesionside = 1 = left
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
                
                ref = mean(ts{1}'); 
                % Voxel-wises correlation 
                zr = fisherz(corr(ref',ts{2}));
                ZR(N.linvoxidxs)=zr; % Expand out to original volume 
                I = N.location2linvoxindxs; 
                i = find(isnan(N.location2linvoxindxs)); 
                I(i)=1; 
                MAP_ZR=ZR(I); 
                MAP_ZR(i)=nan; 
                Si.zr = nanmean(MAP_ZR'); 
                S           = addstruct(S,Si);
            end;
        end;
        save(fullfile(ppDir,'rs_preprocess_voxelwise.mat'),'-struct','S');
        varargout = {S};
        
    case 'PP_voxelwise_analysis' % plotting and statistical analysis of the correlation maps 
        atlasSurfDir = [rootDir '/fmri/surface_caret/fsaverage_sym/LeftHem'];
        Coord = (fullfile(atlasSurfDir,'lh.FLAT.coord'));
        Topo = (fullfile(atlasSurfDir,'lh.CUT.topo'));
        Shape = caret_load(fullfile(atlasSurfDir,'lh.surface_shape'));
        Paint = caret_load(fullfile(atlasSurfDir,'ROI.paint')); 
        B = caret_load('CS.border');
        N   = load(fullfile(roiDir,sprintf('%s_regions_voxelwise.mat','UZP_1005'))); % load subject region definition
        D=load(fullfile(ppDir,'rs_preprocess_voxelwise.mat'));
        figure(1); 
        set(gcf,'PaperPosition',[2 2 12 8]); 
        subplot(3,5,1); 
        [M,d,p]=caret_plotflatmap('coord',Coord,'topo',Topo,'data',Shape.data(:,2),'xlims',[-25 18],'ylims',[-27 42]);
        set(gca,'XTick',[],'YTick',[]); 
        axis equal 
        title('sulcal depth'); 
        subplot(3,5,2); 
        caret_plotflatmap('topo',Topo,'M',M,'data',Paint.data(:,1),'xlims',[-25 18],'ylims',[-27 42]);
        set(gca,'XTick',[],'YTick',[]); 
        axis equal 
        title('ROI'); 
        Data = zeros(163842,size(D.zr,1)); 
        Data(N.location,:) = D.zr'; 
        
        weeks = [1,4,12,24,52]; 
        for i =1:5 
            subplot(3,5,5+i); 
            if (i==1) 
                ylabel('Controls'); 
            end;
            title(sprintf('Week %d', weeks(i))); 
                
            caret_plotflatmap('topo',Topo,'M',M,'data',nanmean(Data(:,D.control & D.week==weeks(i)),2),'cscale',[-0.1 0.5],'border',B.Border); 
            set(gca,'XTick',[],'YTick',[]); 
            axis equal 
            subplot(3,5,10+i); 
            if (i==1) 
                ylabel('patients'); 
            end; 
            caret_plotflatmap('topo',Topo,'M',M,'data',nanmean(Data(:,~D.control & D.week==weeks(i)),2),'cscale',[-0.1 0.5],'border',B.Border); 
            set(gca,'XTick',[],'YTick',[]); 
            axis equal 
        end 
    case 'PP_plotmap'
        M=varargin{1};
        Topo = varargin{2}; 
        
    case 'PP_removePair'                      % some participants have bad correlation pairs, get rid of them
        D = varargin{1};
        
        idx         = (D.r==0 | D.r==1 | D.r==-1);  % edge cases for bad correlation pairs
        D.r(idx)    = nan;
        D.fz_r      = fisherz(D.r);
        
        varargout = {D};
    case 'PP_studySubset'                     % get controls and patients (subcorticals only, more than 1 measurement session)
        D    = load(fullfile(ppDir,'rs_preprocess_test.mat'));
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
        
    case 'NC_getInterHemPattern'              % get pattern of correlations for pairs between hemispheres
        % 0. Load data
        D       = rs_imana('PP_studySubset');
        D       = rs_imana('PP_removePair',D);
        
        % 1. Obtain required data subset
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
                
                % flip hemispheres (left is non-lesioned, right is lesioned)
                if mean(Dw.lesionSide)==1
                    M = M';     % flip to make correlation between pairs for non-lesioned first
                end;
                
                % add data to new structure
                Si      = getrow(Dw,1);
                Si      = rmfield(Si,{'r'});
                Si.M    = nonsym_squareform(M);
                S       = addstruct(S,Si);
            end;
        end;
        save(fullfile(ppDir,'nc_interhem.mat'),'-struct','S');
    case 'NC_getIntraLesionedHemPattern'      % get pattern of correlations for pairs within lesioned hemisphere
        % 0. Load data
        D       = rs_imana('PP_studySubset');
        D       = rs_imana('PP_removePair',D);
        
        % 1. Obtain required data subset
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
                Si      = rmfield(Si,{'r'});
                Si.M    = squareform(M);
                S       = addstruct(S,Si);
            end;
        end;
        save(fullfile(ppDir,'nc_intrahem_lesioned.mat'),'-struct','S');
    case 'NC_getIntraNonlesionedHemPattern'   % get pattern of correlations for pairs within lesioned hemisphere
        % 0. Load data
        D       = rs_imana('PP_studySubset');
        D       = rs_imana('PP_removePair',D);
        
        % 1. Obtain required data subset
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
                Si      = rmfield(Si,{'r'});
                Si.M    = squareform(M);
                S       = addstruct(S,Si);
            end;
        end;
        save(fullfile(ppDir,'nc_intrahem_nonlesioned.mat'),'-struct','S');
    case 'NC_getThemALL'                      % get pattern of correlations for everything
        
        D3 = load(fullfile(ppDir,'nc_intrahem_lesioned.mat'));
        D2 = load(fullfile(ppDir,'nc_intrahem_nonlesioned.mat'));
        D1 = load(fullfile(ppDir,'nc_interhem.mat'));
        
        S = D1;
        S.M =[D1.M, D2.M, D3.M];
        save(fullfile(ppDir,'nc_getThemall.mat'),'-struct','S');
        
    case 'NC_getInterHemPatternSH'            % get pattern of correlations for pairs between hemispheres (split half)
        % 0. Load data
        D       = rs_imana('PP_studySubset');
        D       = rs_imana('PP_removePair',D);
        
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
    case 'NC_getIntraHemLesionedSH'           % get pattern of correlations for pairs between hemispheres (split half)
        % 0. Load data
        D       = rs_imana('PP_studySubset');
        D       = rs_imana('PP_removePair',D);
        
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
    case 'NC_getIntraHemNonlesionedSH'        % get pattern of correlations for pairs between hemispheres (split half)
        % 0. Load data
        D       = rs_imana('PP_studySubset');
        D       = rs_imana('PP_removePair',D);
        
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
        rs_imana('FIG_inter_splithalf_reliability');
        ylim([0 1]);
        ylabel('Fisher Z')
        xlabel('weeks')
        title('Interhemispheric')
        subplot(222);
        rs_imana('FIG_intraLesioned_splithalf_reliability');
        ylim([0 1]);
        ylabel('Fisher Z')
        xlabel('weeks')
        title('Intrahemispheric Les')
        subplot(223);
        rs_imana('FIG_intraNonlesioned_splithalf_reliability');
        ylim([0 1]);
        ylabel('Fisher Z')
        xlabel('weeks')
        title('Intrahemispheric NonLes')
        subplot(224);
        rs_imana('FIG_splithalf_all');
        ylim([0 1]);
        ylabel('Fisher Z')
        xlabel('weeks')
        title('All connections')
        
        keyboard;
        
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
            
            fname   = strcat(pre{1},'_',D.Centre{i},'_',num2str(D.ID(i)),'_',D.Week{i},'.nii');
            sn      = strcat(D.Centre{i},'_',num2str(D.ID(i)));
            fprintf('%d.\t%s\n',D.SN(i),fname);
            
            regFile = fullfile(roiDir,sprintf('%s_regions_all.mat',sn));
            if ~exist(regFile)
                fprintf('No region definition for this subject...\n');
            else
                R   = load(fullfile(roiDir,sprintf('%s_regions_all.mat',sn))); % load subject region definition
                R.R = R.R(inclROI);
                %fDir= fullfile(rsDir,D.subj_name{i},D.Week{i});
                fDir = compareDir;
                v   = spm_vol(fullfile(fDir,fname));                            % load volumes from rs data
                ts  = region_getdata(v,R.R);                                    % get time series for all ROIs
                Bi  = getrow(B,strcmp(B.subj_name,sn) & B.week==D.week(i));
                SSi = getrow(SS,strcmp(SS.subj_name,sn));                       % get gray matter damage from subject_list file
                Li  = getrow(L,strcmp(L.subj_name,sn));                         % get cst damage from DTI
                
                clear Si
                Si.subj_name    = {sn};
                Si.week         = D.week(i);
                Si.lesionType   = find(strcmp(D.LesionLocation(i),lesionType));
                if isempty(find(strcmp(D.LesionLocation(i),lesionType)))
                    Si.lesionType = nan;
                end;
                Si.lesionSide   = D.lesionside(i);
                Si.control      = ~isempty(strfind(sn,'P'));
                if ~isempty(Bi.subj_name)
                    switch(Si.control)
                        case 0
                            Bj = getrow(Bi,Bi.handLabel==4);                    % paretic hand
                        case 1
                            Bj = getrow(Bi,Bi.handLabel==1);                    % non-dom hand
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
                
                % % removing all voxels with no signal
                % ts1 = ts1(:,all(ts1,1));
                % ts2 = ts2(:,all(ts2,1));
                
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
                
                % estimating pair-wise correlations
                r 		= corr(refTS);
                idx 	= find(triu(r,1));
                refFZ 	= fisherz(r(idx));
                refFZ   = refFZ(refFZ>-3 & refFZ<3);
                
                r       = corr(refTS,compTS);
                compFZ  = fisherz(r(:));
                compFZ  = compFZ(compFZ>-3 & compFZ<3);
                
                %                 % estimate pair-wise correlations
                %                 %   - fisherz transform (distribution exists around -3 and 3)
                %                 %   - remove inf
                %                 r       = corr(refTS);
                %                 idx     = any(~isnan(r),1);
                %                 refFZ   = fisherz(1-squareform(1-r(idx,idx)));
                % %                 refFZ   = refFZ(refFZ>-3 & refFZ<3);
                
                %                 r       = corr(refTS,compTS);
                %                 compFZ  = fisherz(r(~isnan(r)));
                %                 compFZ  = compFZ(:);
                %                 % compFZ  = compFZ(compFZ>-3 & compFZ<3);
                
                Si.compFZ   = nanmean(compFZ);
                Si.refFZ    = nanmean(refFZ);
                
                Si.sumnan_compFZ   = sum(isnan(compFZ));
                Si.sumnan_refFZ    = sum(isnan(refFZ));
                
                fprintf('comp: %d, ref %d\n',Si.sumnan_compFZ,Si.sumnan_refFZ);
                %                 Si.refCon   = mean(Si.compFZ)./mean(Si.refFZ);
                
                S           = addstruct(S,Si);
            end;
        end;
        varargout = {S};
        save(fullfile(ppDir,'golestani2013_preprocess_new.mat'),'-struct','S');
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
                %fDir= fullfile(rsDir,D.subj_name{i},D.Week{i});
                fDir = compareDir;
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
                Si.lesionSide   = D.lesionside(i);                          %lesionside = 1
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
                % cComp(isinf(cComp)) = nan;
                cComp = fisherinv(nanmean(cComp,1))';
                
                % estimate 95% percentile of correlations
                
                allcorr_nan = [cRef;cComp];
                allcorr = allcorr_nan(~isnan(allcorr_nan));
                pr      = prctile(allcorr,95);
                
                if sum(isinf([cRef;cComp]))>0
                    keyboard;
                end;
                
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
        save(fullfile(ppDir,'park2011_preprocess_new.mat'),'-struct','S');
        
        
        %% Analysis for the Paper
        % New Analysis with Euclidian distances
        % Aufbau: last case with Vorsilbe FIG is always the figure combing the preceeding cases
        % first case has always more explanatory comments
    case 'FIG_M1'
        D       = rs_imana('PP_studySubset');
        D       = rs_imana('PP_removePair',D);
        
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
        
        keyboard;
    case 'FIG_Park'
        D   = load(fullfile(ppDir,'park2011_preprocess_new.mat'));
        
        
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
    case 'FIG_RelCon'
        D   = load(fullfile(ppDir,'golestani2013_preprocess_new.mat'));
        
        % 1. get controls
        Dc      = getrow(D,D.control==1);
        [x,sn]  = pivottable(Dc.subj_name,[],Dc.week,'length(unique(x))');
        Dc      = getrow(Dc,ismember(Dc.subj_name,sn(x>1)));
        
        Dp      = getrow(D,D.lesionType==3);
        [x,sn]  = pivottable(Dp.subj_name,[],Dp.week,'length(unique(x))');
        Dp      = getrow(Dp,ismember(Dp.subj_name,sn(x>1)));
        
        T = addstruct(Dc,Dp);
        varargout = {T};
        
        T.RelCon = T.compFZ./T.refFZ;
        
        h = make_figure;
        s   = style_sheet(sty_2grp_ac,'leg',{'patient','control'},'leglocation','northoutside','errorcolor',sty_2grp_ac);
        
        title('RelCon for interhemispheric Sensorimotor cortex (S1M1-S1M1)')
        lineplot(T.week,T.RelCon,'plotfcn','nanmean','errorfcn','stderr','split',T.control,'CAT',s(1).CAT,s(1).PLOT{:});
        ylabel('RelCon')
        xlabel('weeks')
        
        % correlation between RelCon and FM
        % needs to be done for hand strength and ARAT too
        
        %      Tp = getrow(T,T.control==0);
        %
        %      T1     = getrow(Tp,Tp.week==1);
        %      T4     = getrow(Tp,Tp.week==4);
        %      T12    = getrow(Tp,Tp.week==12);
        %      T24    = getrow(Tp,Tp.week==24);
        %      T52    = getrow(Tp,Tp.week==52);
        %
        %      [r1,p1] = corr(T1.refFZ,T1.FM,'row','complete');
        %      [r2,p2] = corr(T4.refFZ,T4.FM,'row','complete');
        %      [r3,p3] = corr(T12.refFZ,T12.FM,'row','complete');
        %      [r4,p4] = corr(T24.refFZ,T24.FM,'row','complete');
        %      [r5,p5] = corr(T52.refFZ,T52.FM,'row','complete');
        
        DD =     tapply(T,{'subj_name','week','control'},{'RelCon'});
        dsave(fullfile(statsDir,'RelCon_stat.dat'),DD);
        
        keyboard;
        
    case 'FIG_interhemSomat'                       % somatotopic ordering of correlations between hemispheres
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
    case 'control_patient_test' % Implements Patient vs. control Null-hypothesis and Equivalence test
        % Examples:
        % rs_imana('control_patient_test','numIter',10000,'euc_true',[0 0.37],'connect','lesioned','week',1);
        % rs_imana('control_patient_test','numIter',10000,'euc_true',[0 1.4],'connect','all','week',1);
        
        week = 1;
        numIter = 1000;
        euc_true = [0 1.4]; % Which true Eucledian distance to test
        connect = 'lesioned';
        vararginoptions(varargin,{'week','numIter','euc_true','connect'});
        
        % determine the correct file
        switch(connect)
            case 'lesioned'
                file = 'nc_intrahem_lesioned.mat';
            case 'nonlesioned'
                file = 'nc_intrahem_nonlesioned.mat';
            case 'interhem'
                file = 'nc_interhem.mat';
            case 'all'
                file = 'nc_getThemall.mat';
        end;
        
        
        D = load(fullfile(ppDir,file)); % _joern
        D  = getrow(D,D.week==week);
        D.res = bsxfun(@minus,D.M,mean(D.M)); % Get the residuals
        numConnect = size(D.M,2);  % Number of connections
        
        numC = sum(D.control);
        numP = sum(~D.control);
        N = numC + numP;
        
        % Calculate the real difference
        diff = sum(D.res(D.control,:))/numC-sum(D.res(~D.control,:))/numP;
        euc_emp = sqrt(sum(diff.^2));
        
        % Estimate a approximate effect size for indepdent connections
        res(D.control,:) = bsxfun(@minus,D.res(D.control,:),mean(D.res(D.control,:))); % Get the residuals
        res(~D.control,:) = bsxfun(@minus,D.res(~D.control,:),mean(D.res(~D.control,:))); % Get the residuals
        resms = std(res);
        effect_size_emp = euc_emp/sqrt(numConnect)/mean(resms); % univariate effect size
        
        PP=[];
        for k=1:length(euc_true)
            P=[];
            for j=1:numIter
                % Make a new assignment of controls and patients
                contT = false(N,1);
                a=sample_wor([1:N],numC);
                contT(a)=true;
                
                % Add a real effect to the patients
                realDiff = normrnd(0,1,1,numConnect);
                realDiff=realDiff/sqrt(sum(realDiff.^2))*euc_true(k);
                Data = D.res;
                Data(~contT,:) = bsxfun(@plus,Data(~contT,:),realDiff);
                
                % Now dtermine the empirically measured distance under this
                % true effect (larger than true effect)
                diff = sum(Data(contT,:))/numC-sum(Data(~contT,:))/numP;
                P.euc_sim(j,1) = sqrt(sum(diff.^2));
                
                % evalauate and record the results
                P.euc_true(j,1) = euc_true(k);
            end;
            PP=addstruct(PP,P);
            EUC{k}=P.euc_sim;
            if (euc_true(k)==0)             % Evaluate against Null hypothesis
                fprintf('%s Week: %d\n',connect,week);
                fprintf('delta_pattern = %2.3f; Empiricial univariate effect = %2.3f\n', euc_emp,effect_size_emp);
                p_val = sum(P.euc_sim>=euc_emp)/numIter;
                % Critical value for 5% significance (for gray interval)
                crit_val =  prctile(P.euc_sim,95);
                fprintf('Against Null: p-value = %2.3f Crtical-val = %2.3f\n', p_val,crit_val);
            else                            % Evaluate against Alternative hypothesis
                p_val = sum(P.euc_sim<=euc_emp)/numIter;
                effect = euc_true(k)/sqrt(numConnect)/mean(resms); % univariate effect size
                fprintf('Against Alternative: delta: %2.3f, univariate effect: %2.3f, p-value= %2.3f\n',euc_true(k),effect,p_val);
            end;
        end;
        nhist(EUC,'noerror','smooth','color','sequential','samebins');
        drawline(euc_emp,'dir','vert');
        drawline(euc_true(1),'dir','vert','color','b');
        drawline(euc_true(2),'dir','vert','color','r');
        set_graphics(gcf,'ylabel',{'frequency'}, 'xlabel',{'Euclidian distance'},'fontsize_all',14);
        set(gcf,'PaperPosition',[2 2 6 4]);
        wysiwyg;
        varargout={PP};
    case 'acute_control_patient_effectsize'
        numRand = 1000; % Randomization for permutation test
        numIter = 2000;   % Number iterations for each effect size
        D = load(fullfile(ppDir,'nc_getThemall.mat')); % _joern
        D  = getrow(D,D.week==1);                                          % Dw needs to be changed manually if you want to look at cross-sectional comparison of other weeks
        
        D.res = bsxfun(@minus,D.M,mean(D.M)); % Get the residuals
        
        numC = sum(D.control);
        numP = sum(~D.control);
        N = numC + numP;
        
        PP=[];
        euc_true = [0:0.1:3];
        for k=1:length(euc_true)
            P=[];
            for j=1:numIter
                % Make a new assignment of controls and patients
                contT = false(N,1);
                a=sample_wor([1:N],numC);
                contT(a)=true;
                
                % Add a real effect to the patients
                realDiff = normrnd(0,1,1,45);
                realDiff=realDiff/sqrt(sum(realDiff.^2))*euc_true(k);
                Data = D.res;
                Data(~contT,:) = bsxfun(@plus,Data(~contT,:),realDiff);
                
                % Now dtermine the empirically measured distance under this
                % true effect (larger than true effect)
                diff = sum(Data(contT,:))/numC-sum(Data(~contT,:))/numP;
                euc_real = sqrt(sum(diff.^2));
                euc = zeros(numRand,1);
                
                % Now run boostrap test
                for i=1:numRand
                    contR = false(N,1);
                    a=sample_wor([1:N],numC);
                    contR(a)=true;
                    diff = sum(Data(contR,:))/numC-sum(Data(~contR,:))/numP;
                    euc(i) = sqrt(sum(diff.^2));
                end;
                % evalauate and record the results
                P.p(j,1) = sum(euc>euc_real)/numRand;
                P.euc_real(j,1) = euc_real;
                P.euc_true(j,1) = euc_true(k);
            end;
            PP.euc_true(k,1) = euc_true(k);
            PP.euc_mean(k,1) = mean(P.euc_real);
            PP.euc_std(k,1) = std(P.euc_real);
            PP.sig5(k,1) = sum(P.p<0.05)/numIter;
            PP.sig1(k,1) = sum(P.p<0.01)/numIter;
        end;
        varargout={PP};
    case 'acute_control_patient_effectsize_plot'
        D=varargin{1};
        plot(D.euc_true,D.sig5,'r','LineWidth',2);
        hold on;
        plot(D.euc_true,D.sig1,'b','LineWidth',2);
        xlabel('TrueDistance')
        ylabel('Proportion positive');
        drawline([0.01 0.05 0.8],'dir','horz');
        set(gca,'Box','off');
    case 'acute_control_patient_all_connections'   % acute crossectional comparison patients versus controls
        
        D = load(fullfile(ppDir,'nc_getThemall.mat')); % _joern
        
        Dw  = getrow(D,D.week==1);                                          % Dw needs to be changed manually if you want to look at cross-sectional comparison of other weeks
        %Dw  = getrow(Dw,Dw.FM>40); %when only looking at mild patients
        
        Dc  = getrow(Dw,Dw.control==1);
        Dp  = getrow(Dw,Dw.control==0);
        %Dp  = getrow(Dp,Dp.FM<=40); %when only looking at severe patients
        
        
        euc_real = pdist2((mean(Dc.M)), (mean(Dp.M)),'euclidean');
        
        
        % permutation
        P = [];
        h = waitbar(0,'Please wait, Bootstraping the Data');
        NumPer=10000;
        
        %CAVE for severe patients
        Dcomb = [Dp,Dc];
        combo = struct;  %final structure
        for field = fieldnames(Dcomb)'
            fname = field{1};
            combo.(fname) = vertcat(Dcomb.(fname));
        end
        
        
        %[~,~,Dw.SN] = unique(Dw.subj_name);
        [~,~,combo.SN] = unique(combo.subj_name);
        
        for j=1:NumPer
            
            fake_control   = datasample(unique(combo.SN),11);
            fake_subject   = combo.SN(~ismember(combo.SN,fake_control));
            
            F_acute_control = combo.M(fake_control,:);
            F_acute_patient = combo.M(fake_subject,:);
            
            Pi.euc_fake              = pdist2(mean(F_acute_control), mean(F_acute_patient),'euclidean');
            
            Pi.j = j;
            
            h=waitbar(j / NumPer);
            P = addstruct(P,Pi);
            
        end
        
        Numb_exceed       = (euc_real<P.euc_fake);
        EUC_per_all       = (sum(Numb_exceed))/10000;
        keyboard
        
        nhist(P.euc_fake,'noerror','smooth','color','sequential','samebins');
        drawline(euc_real,'dir','vert','linestyle','--','color',[1 0 0])
        % 5% & 95% Intervalle
        y = prctile(P.euc_fake, [2.5 97.5]);
        drawline(y,'dir','vert','linestyle',':','color',[0 0 0])
        set_graphics(gcf,'ylabel',{'frequency'}, 'xlabel',{'Euclidian distance'},'fontsize_all',14);
        %     dim = [.5 .5 .3 .3];
        %     str = EUC_per_all;
        %     annotation('textbox',dim,'String',str,'FitBoxToText','on');
    case 'acute_control_patient_inter'
        D = load(fullfile(ppDir,'nc_interhem.mat'));
        
        Dw  = getrow(D,D.week==1);
        %Dw  = getrow(Dw,Dw.FM>40);
        
        Dc  = getrow(Dw,Dw.control==1);
        Dp  = getrow(Dw,Dw.control==0);
        Dp  = getrow(Dp,Dp.FM<=40); %when only looking at severe patients
        
        D.euc_real = pdist2((mean(Dc.M)), (mean(Dp.M)),'euclidean');
        
        % permutation
        P = [];
        h = waitbar(0,'Please wait, Bootstraping the Data');
        NumPer=10000;
        
        %[~,~,Dw.SN] = unique(Dw.subj_name);
        %CAVE for severe patients
        Dcomb = [Dp,Dc];
        combo = struct;  %final structure
        for field = fieldnames(Dcomb)'
            fname = field{1};
            combo.(fname) = vertcat(Dcomb.(fname));
        end
        [~,~,combo.SN] = unique(combo.subj_name);
        
        for j=1:NumPer
            
            fake_control   = datasample(unique(combo.SN),11);
            fake_subject   = combo.SN(~ismember(combo.SN,fake_control));
            
            F_acute_control = combo.M(fake_control,:);
            F_acute_patient = combo.M(fake_subject,:);
            
            Pi.euc_fake              = pdist2(mean(F_acute_control), mean(F_acute_patient),'euclidean');
            
            Pi.j = j;
            
            h=waitbar(j / NumPer);
            P = addstruct(P,Pi);
            
        end
        
        Numb_exceed       = (D.euc_real<P.euc_fake);
        EUC_per_inter       = (sum(Numb_exceed))/10000;
        keyboard
        
        nhist(P.euc_fake,'noerror','smooth','color','sequential','samebins');
        drawline(D.euc_real,'dir','vert','linestyle','--','color',[1 0 0])
        % 5% & 95% Intervalle
        y = prctile(P.euc_fake, [2.5 97.5]);
        drawline(y,'dir','vert','linestyle',':','color',[0 0 0])
        set_graphics(gcf,'ylabel',{'frequency'}, 'xlabel',{'mean distance'},'fontsize_all',14);
        %     dim = [.8 .6 .3 .3];
        %     str = EUC_per_inter;
        %     annotation('textbox',dim,'String',str,'FitBoxToText','on');
    case 'acute_control_patient_intra_les'
        D = load(fullfile(ppDir,'nc_intrahem_lesioned.mat'));
        
        Dw  = getrow(D,D.week==1);
        %Dw  = getrow(Dw,Dw.FM>40);
        
        Dc  = getrow(Dw,Dw.control==1);
        Dp  = getrow(Dw,Dw.control==0);
        Dp  = getrow(Dp,Dp.FM<=40); %when only looking at severe patients
        
        D.euc_real = pdist2((mean(Dc.M)), (mean(Dp.M)),'euclidean');
        
        % permutation
        P = [];
        h = waitbar(0,'Please wait, Bootstraping the Data');
        NumPer=10000;
        
        %[~,~,Dw.SN] = unique(Dw.subj_name);
        %CAVE for severe patients
        Dcomb = [Dp,Dc];
        combo = struct;  %final structure
        for field = fieldnames(Dcomb)'
            fname = field{1};
            combo.(fname) = vertcat(Dcomb.(fname));
        end
        
        [~,~,combo.SN] = unique(combo.subj_name);
        
        for j=1:NumPer
            
            fake_control   = datasample(unique(combo.SN),11);
            fake_subject   = combo.SN(~ismember(combo.SN,fake_control));
            
            F_acute_control = combo.M(fake_control,:);
            F_acute_patient = combo.M(fake_subject,:);
            
            Pi.euc_fake              = pdist2(mean(F_acute_control), mean(F_acute_patient),'euclidean');
            
            Pi.j = j;
            
            h=waitbar(j / NumPer);
            P = addstruct(P,Pi);
            
        end
        
        Numb_exceed             = (D.euc_real<P.euc_fake);
        EUC_per_intra_les       = (sum(Numb_exceed))/10000;
        keyboard
        
        nhist(P.euc_fake,'noerror','smooth','color','sequential','samebins');
        drawline(D.euc_real,'dir','vert','linestyle','--','color',[1 0 0])
        % 5% & 95% Intervalle
        y = prctile(P.euc_fake, [2.5 97.5]);
        drawline(y,'dir','vert','linestyle',':','color',[0 0 0])
        set_graphics(gcf,'ylabel',{'frequency'}, 'xlabel',{'mean distance'},'fontsize_all',14);
        %     dim = [.5 .5 .3 .3];
        %     str = EUC_per_intra_les;
        %     annotation('textbox',dim,'String',str,'FitBoxToText','on');
    case 'acute_control_patient_intra_non'
        
        D = load(fullfile(ppDir,'nc_intrahem_nonlesioned.mat'));
        
        Dw  = getrow(D,D.week==1);
        %Dw  = getrow(Dw,Dw.FM>40);
        
        Dc  = getrow(Dw,Dw.control==1);
        Dp  = getrow(Dw,Dw.control==0);
        Dp  = getrow(Dp,Dp.FM<=40); %when only looking at severe patients
        
        D.euc_real = pdist2((mean(Dc.M)), (mean(Dp.M)),'euclidean');
        
        % permutation
        P = [];
        h = waitbar(0,'Please wait, Bootstraping the Data');
        NumPer=10000;
        
        %[~,~,Dw.SN] = unique(Dw.subj_name);
        %CAVE for severe patients
        Dcomb = [Dp,Dc];
        combo = struct;  %final structure
        for field = fieldnames(Dcomb)'
            fname = field{1};
            combo.(fname) = vertcat(Dcomb.(fname));
        end
        
        [~,~,combo.SN] = unique(combo.subj_name);
        
        for j=1:NumPer
            
            fake_control   = datasample(unique(combo.SN),11);
            fake_subject   = combo.SN(~ismember(combo.SN,fake_control));
            
            F_acute_control = combo.M(fake_control,:);
            F_acute_patient = combo.M(fake_subject,:);
            
            Pi.euc_fake              = pdist2(mean(F_acute_control), mean(F_acute_patient),'euclidean');
            
            Pi.j = j;
            
            h=waitbar(j / NumPer);
            P = addstruct(P,Pi);
            
        end
        
        Numb_exceed             = (D.euc_real<P.euc_fake);
        EUC_per_intra_non       = (sum(Numb_exceed))/10000;
        keyboard
        
        nhist(P.euc_fake,'noerror','smooth','color','sequential','samebins');
        drawline(D.euc_real,'dir','vert','linestyle','--','color',[1 0 0])
        % 5% & 95% Intervalle
        y = prctile(P.euc_fake, [2.5 97.5]);
        drawline(y,'dir','vert','linestyle',':','color',[0 0 0])
        set_graphics(gcf,'ylabel',{'frequency'}, 'xlabel',{'mean distance'},'fontsize_all',14);
        %     dim = [.5 .5 .3 .3];
        %     str = EUC_per_intra_non;
        %     annotation('textbox',dim,'String',str,'FitBoxToText','on');
    case 'FIG_acute_control_patient'
        subplot(3,3,[7 2])
        rs_imana('acute_control_patient_all_connections');
        title ('full connectivity pattern')
        
        subplot(3,3,3)
        rs_imana('acute_control_patient_inter');
        title ('interhemispheric connectivity pattern')
        
        subplot(3,3,6)
        rs_imana('acute_control_patient_intra_les');
        title ('intrahemispheric lesioned connectivity pattern')
        
        subplot(3,3,9)
        rs_imana('acute_control_patient_intra_non');
        title ('intrahemispheric nonlesioned connectivity pattern')
        keyboard;
        
    case 'single_connections'
        
        D = load(fullfile(ppDir,'nc_getThemall_jorn.mat'));
        
        Dw  = getrow(D,D.week==52);                                          % Dw needs to be changed manually if you want to look at cross-sectional comparison of other weeks
        %Dw  = getrow(Dw,Dw.FM>40); %when only looking at mild patients
        
        Dc  = getrow(Dw,Dw.control==1);
        Dp  = getrow(Dw,Dw.control==0);
        
        
        for i = 1:length(D.M(1,:))
            x = pivottable(Dc.subj_name,Dc.week,Dc.M(:,i),'nanmean');
            y = pivottable(Dp.subj_name,Dp.week,Dp.M(:,i),'nanmean');
            
            ttest(x(:,1),y(:,1), 2,'independent');
            
        end;
        
    case 'acute_variability_all'                   % acute variability comparison patients versus controls
        D = load(fullfile(ppDir,'nc_getThemall_jorn.mat'));
        
        Dw  = getrow(D,D.week==52);                                          % Dw needs to be changed manually if you want to look at cross-sectional comparison of other weeks
        
        Dc  = getrow(Dw,Dw.control==1);
        Dp  = getrow(Dw,Dw.control==0);
        
        C = [];
        P = [];
        leaveoutCAvg = nan(length(Dc.subj_name),45);
        leaveoutPAvg = nan(length(Dp.subj_name),45);
        
        %job is to calculate the distance between the average and each
        %subject to get a measure for Variability
        
        idx = 1;
        for c=unique(Dc.subj_name)'
            isC     = getrow(Dc,ismember(Dc.subj_name,c));
            notC    = getrow(Dc,~ismember(Dc.subj_name,c));
            % add control correlations with controls
            Ci.subj_name    = c;
            % calculate leave-one-out average of all controls
            leaveoutCAvg(idx,:)  = nanmean(notC.M,1);
            Ci.euc              = pdist2(isC.M,leaveoutCAvg(idx,:),'euclidean');
            idx             = idx + 1;
            C               = addstruct(C,Ci);
        end;
        
        idx = 1;
        for p=unique(Dp.subj_name)'
            isP     = getrow(Dp,ismember(Dp.subj_name,p));
            notP    = getrow(Dp,~ismember(Dp.subj_name,p));
            % add patients correlations with patients
            Pi.subj_name    = p;
            % calculate leave-one-out average of all patients
            leaveoutPAvg(idx,:)  = nanmean(notP.M,1);
            Pi.euc              = pdist2(isP.M,leaveoutPAvg(idx,:),'euclidean');
            idx             = idx + 1;
            P               = addstruct(P,Pi);
        end;
        
        %
        EUC_real = (nanmean(P.euc))-(nanmean(C.euc));
        
        % Permutation
        Per = [];
        h = waitbar(0,'Please wait, Bootstraping the Data');
        NumPer=10000;
        
        [~,~,Dw.SN] = unique(Dw.subj_name);
        
        FP = [];
        FC = [];
        for j=1:NumPer
            
            fake_control   = datasample(unique(Dw.SN),11);
            fake_subject   = Dw.SN(~ismember(Dw.SN,fake_control));
            %             fake_subject   = datasample(unique(Dw.SN),16);
            
            F_control = Dw.M(fake_control,:);
            F_patient = Dw.M(fake_subject,:);
            
            leaveoutCAvg = nan(length(F_control),45);
            leaveoutPAvg = nan(length(F_patient),45);
            
            fcidx = 1;
            for fc=1:length(fake_control)
                isC         = F_control(fcidx,:);
                notC        = F_control((1:end ~=fcidx),:);
                % calculate leave-one-out average of all controls
                FCi.euc     = pdist2(isC,nanmean(notC),'euclidean');
                fcidx       = fcidx + 1;
                FC          = addstruct(FC,FCi);
            end;
            
            fpidx = 1;
            for fp=1:length(fake_subject)
                isP         = F_patient(fpidx,:);
                notP        = F_patient((1:end ~=fpidx),:);
                % calculate leave-one-out average of all patients
                FPi.euc     = pdist2(isP,nanmean(notP),'euclidean');
                fpidx       = fpidx + 1;
                FP          = addstruct(FP,FPi);
            end;
            
            PP.EUC_fake = (nanmean(FP.euc))-(nanmean(FC.euc));
            
            PP.j = j;
            
            h=waitbar(j / NumPer);
            Per = addstruct(Per,PP);
        end
        Numb_exceed             = (EUC_real<Per.EUC_fake);
        EUC_per_all      = (sum(Numb_exceed))/10000;
        
        nhist(Per.EUC_fake,'noerror','smooth','color','sequential','samebins');
        drawline(EUC_real,'dir','vert','linestyle','--','color',[1 0 0])
        % 5% & 95% Intervalle
        y = prctile(Per.EUC_fake, [2.5 97.5]);
        drawline(y,'dir','vert','linestyle',':','color',[0 0 0])
        set_graphics(gcf,'ylabel',{'frequency'}, 'xlabel',{'mean distance'},'fontsize_all',14);
        xlim([-0.04 0.55])
        keyboard;
    case 'acute_variability_inter'
        D = load(fullfile(ppDir,'nc_interhem.mat'));
        
        Dw  = getrow(D,D.week==4);
        
        Dc  = getrow(Dw,Dw.control==1);
        Dp  = getrow(Dw,Dw.control==0);
        
        C = [];
        P = [];
        leaveoutCAvg = nan(length(Dc.subj_name),25);
        leaveoutPAvg = nan(length(Dp.subj_name),25);
        
        idx = 1;
        for c=unique(Dc.subj_name)'
            isC     = getrow(Dc,ismember(Dc.subj_name,c));
            notC    = getrow(Dc,~ismember(Dc.subj_name,c));
            % add control correlations with controls
            Ci.subj_name    = c;
            % calculate leave-one-out average of all controls
            leaveoutCAvg(idx,:)  = nanmean(notC.M,1);
            Ci.euc              = pdist2(isC.M,leaveoutCAvg(idx,:),'euclidean');
            idx             = idx + 1;
            C               = addstruct(C,Ci);
        end;
        
        idx = 1;
        for p=unique(Dp.subj_name)'
            isP     = getrow(Dp,ismember(Dp.subj_name,p));
            notP    = getrow(Dp,~ismember(Dp.subj_name,p));
            % add patients correlations with patients
            Pi.subj_name    = p;
            % calculate leave-one-out average of all patients
            leaveoutPAvg(idx,:)  = nanmean(notP.M,1);
            Pi.euc              = pdist2(isP.M,leaveoutPAvg(idx,:),'euclidean');
            idx             = idx + 1;
            P               = addstruct(P,Pi);
        end;
        
        EUC_real = (mean(P.euc))-(mean(C.euc));
        
        % Permutation
        Per = [];
        h = waitbar(0,'Please wait, Bootstraping the Data');
        NumPer=10000;
        
        [~,~,Dw.SN] = unique(Dw.subj_name);
        
        FP = [];
        FC = [];
        for j=1:NumPer
            
            fake_control   = datasample(unique(Dw.SN),11);
            fake_subject   = Dw.SN(~ismember(Dw.SN,fake_control));
            %             fake_subject   = datasample(unique(Dw.SN),16);
            
            F_control = Dw.M(fake_control,:);
            F_patient = Dw.M(fake_subject,:);
            
            leaveoutCAvg = nan(length(F_control),45);
            leaveoutPAvg = nan(length(F_patient),45);
            
            fcidx = 1;
            for fc=1:length(fake_control)
                isC     = F_control(fcidx,:);
                notC     = F_control((1:end ~=fcidx),:);
                % calculate leave-one-out average of all controls
                FCi.euc               = pdist2(isC,nanmean(notC),'euclidean');
                fcidx             = fcidx + 1;
                FC               = addstruct(FC,FCi);
            end;
            
            fpidx = 1;
            for fp=1:length(fake_subject)
                isP         = F_patient(fpidx,:);
                notP        = F_patient((1:end ~=fpidx),:);
                % calculate leave-one-out average of all patients
                FPi.euc     = pdist2(isP,nanmean(notP),'euclidean');
                fpidx       = fpidx + 1;
                FP          = addstruct(FP,FPi);
            end;
            
            PP.EUC_fake = (mean(FP.euc))-(mean(FC.euc));
            
            PP.j = j;
            
            h=waitbar(j / NumPer);
            Per = addstruct(Per,PP);
        end
        
        Numb_exceed             = (EUC_real<Per.EUC_fake);
        EUC_per_inter       = (sum(Numb_exceed))/10000;
        
        nhist(Per.EUC_fake,'noerror','smooth','color','sequential','samebins');
        drawline(EUC_real,'dir','vert','linestyle','--','color',[1 0 0])
        y = prctile(Per.EUC_fake, [2.5 97.5]);
        drawline(y,'dir','vert','linestyle',':','color',[0 0 0])
        set_graphics(gcf,'ylabel',{'frequency'}, 'xlabel',{'mean distance'},'fontsize_all',14);
        xlim([-0.04 0.55])
        keyboard;
    case 'acute_variability_intra_les'
        D = load(fullfile(ppDir,'nc_intrahem_lesioned.mat'));
        Dw  = getrow(D,D.week==4);
        
        Dc  = getrow(Dw,Dw.control==1);
        Dp  = getrow(Dw,Dw.control==0);
        
        C = [];
        P = [];
        leaveoutCAvg = nan(length(Dc.subj_name),10);
        leaveoutPAvg = nan(length(Dp.subj_name),10);
        
        idx = 1;
        for c=unique(Dc.subj_name)'
            isC     = getrow(Dc,ismember(Dc.subj_name,c));
            notC    = getrow(Dc,~ismember(Dc.subj_name,c));
            % add control correlations with controls
            Ci.subj_name    = c;
            % calculate leave-one-out average of all controls
            leaveoutCAvg(idx,:)  = nanmean(notC.M,1);
            Ci.euc              = pdist2(isC.M,leaveoutCAvg(idx,:),'euclidean');
            idx             = idx + 1;
            C               = addstruct(C,Ci);
        end;
        
        idx = 1;
        for p=unique(Dp.subj_name)'
            isP     = getrow(Dp,ismember(Dp.subj_name,p));
            notP    = getrow(Dp,~ismember(Dp.subj_name,p));
            % add patients correlations with patients
            Pi.subj_name    = p;
            % calculate leave-one-out average of all patients
            leaveoutPAvg(idx,:)  = nanmean(notP.M,1);
            Pi.euc              = pdist2(isP.M,leaveoutPAvg(idx,:),'euclidean');
            idx             = idx + 1;
            P               = addstruct(P,Pi);
        end;
        
        EUC_real = (mean(P.euc))-(mean(C.euc));
        
        % Permutation
        Per = [];
        h = waitbar(0,'Please wait, Bootstraping the Data');
        NumPer=10000;
        
        [~,~,Dw.SN] = unique(Dw.subj_name);
        
        FP = [];
        FC = [];
        for j=1:NumPer
            
            fake_control   = datasample(unique(Dw.SN),11);
            fake_subject   = Dw.SN(~ismember(Dw.SN,fake_control));
            %             fake_subject   = datasample(unique(Dw.SN),16);
            
            F_control = Dw.M(fake_control,:);
            F_patient = Dw.M(fake_subject,:);
            
            leaveoutCAvg = nan(length(F_control),45);
            leaveoutPAvg = nan(length(F_patient),45);
            
            fcidx = 1;
            for fc=1:length(fake_control)
                isC     = F_control(fcidx,:);
                notC     = F_control((1:end ~=fcidx),:);
                % calculate leave-one-out average of all controls
                FCi.euc               = pdist2(isC,nanmean(notC),'euclidean');
                fcidx             = fcidx + 1;
                FC               = addstruct(FC,FCi);
            end;
            
            fpidx = 1;
            for fp=1:length(fake_subject)
                isP         = F_patient(fpidx,:);
                notP        = F_patient((1:end ~=fpidx),:);
                % calculate leave-one-out average of all patients
                FPi.euc     = pdist2(isP,nanmean(notP),'euclidean');
                fpidx       = fpidx + 1;
                FP          = addstruct(FP,FPi);
            end;
            
            PP.EUC_fake = (mean(FP.euc))-(mean(FC.euc));
            
            PP.j = j;
            
            h=waitbar(j / NumPer);
            Per = addstruct(Per,PP);
        end
        Numb_exceed             = (EUC_real<Per.EUC_fake);
        EUC_per_intra_les       = (sum(Numb_exceed))/10000;
        
        nhist(Per.EUC_fake,'noerror','smooth','color','sequential','samebins');
        drawline(EUC_real,'dir','vert','linestyle','--','color',[1 0 0])
        y = prctile(Per.EUC_fake, [2.5 97.5]);
        drawline(y,'dir','vert','linestyle',':','color',[0 0 0])
        set_graphics(gcf,'ylabel',{'frequency'}, 'xlabel',{'mean distance'},'fontsize_all',14);
        xlim([-0.04 0.55])
        keyboard;
    case 'acute_variability_intra_non'
        D = load(fullfile(ppDir,'nc_intrahem_nonlesioned.mat'));
        
        Dw  = getrow(D,D.week==4);
        
        Dc  = getrow(Dw,Dw.control==1);
        Dp  = getrow(Dw,Dw.control==0);
        
        C = [];
        P = [];
        leaveoutCAvg = nan(length(Dc.subj_name),10);
        leaveoutPAvg = nan(length(Dp.subj_name),10);
        
        idx = 1;
        for c=unique(Dc.subj_name)'
            isC     = getrow(Dc,ismember(Dc.subj_name,c));
            notC    = getrow(Dc,~ismember(Dc.subj_name,c));
            % add control correlations with controls
            Ci.subj_name    = c;
            % calculate leave-one-out average of all controls
            leaveoutCAvg(idx,:)  = nanmean(notC.M,1);
            Ci.euc              = pdist2(isC.M,leaveoutCAvg(idx,:),'euclidean');
            idx             = idx + 1;
            C               = addstruct(C,Ci);
        end;
        
        idx = 1;
        for p=unique(Dp.subj_name)'
            isP     = getrow(Dp,ismember(Dp.subj_name,p));
            notP    = getrow(Dp,~ismember(Dp.subj_name,p));
            % add patients correlations with patients
            Pi.subj_name    = p;
            % calculate leave-one-out average of all patients
            leaveoutPAvg(idx,:)  = nanmean(notP.M,1);
            Pi.euc              = pdist2(isP.M,leaveoutPAvg(idx,:),'euclidean');
            idx             = idx + 1;
            P               = addstruct(P,Pi);
        end;
        
        EUC_real = (mean(P.euc))-(mean(C.euc));
        
        % Permutation
        Per = [];
        h = waitbar(0,'Please wait, Bootstraping the Data');
        NumPer=10000;
        
        [~,~,Dw.SN] = unique(Dw.subj_name);
        
        FP = [];
        FC = [];
        for j=1:NumPer
            
            fake_control   = datasample(unique(Dw.SN),11);
            fake_subject   = Dw.SN(~ismember(Dw.SN,fake_control));
            %             fake_subject   = datasample(unique(Dw.SN),16);
            
            F_control = Dw.M(fake_control,:);
            F_patient = Dw.M(fake_subject,:);
            
            leaveoutCAvg = nan(length(F_control),45);
            leaveoutPAvg = nan(length(F_patient),45);
            
            fcidx = 1;
            for fc=1:length(fake_control)
                isC     = F_control(fcidx,:);
                notC     = F_control((1:end ~=fcidx),:);
                % calculate leave-one-out average of all controls
                FCi.euc               = pdist2(isC,nanmean(notC),'euclidean');
                fcidx             = fcidx + 1;
                FC               = addstruct(FC,FCi);
            end;
            
            fpidx = 1;
            for fp=1:length(fake_subject)
                isP         = F_patient(fpidx,:);
                notP        = F_patient((1:end ~=fpidx),:);
                % calculate leave-one-out average of all patients
                FPi.euc     = pdist2(isP,nanmean(notP),'euclidean');
                fpidx       = fpidx + 1;
                FP          = addstruct(FP,FPi);
            end;
            
            PP.EUC_fake = (mean(FP.euc))-(mean(FC.euc));
            
            PP.j = j;
            
            h=waitbar(j / NumPer);
            Per = addstruct(Per,PP);
        end
        Numb_exceed             = (EUC_real<Per.EUC_fake);
        EUC_per_intra_non       = (sum(Numb_exceed))/10000;
        
        nhist(Per.EUC_fake,'noerror','smooth','color','sequential','samebins');
        drawline(EUC_real,'dir','vert','linestyle','--','color',[1 0 0])
        y = prctile(Per.EUC_fake, [2.5 97.5]);
        drawline(y,'dir','vert','linestyle',':','color',[0 0 0])
        set_graphics(gcf,'ylabel',{'frequency'}, 'xlabel',{'mean distance'},'fontsize_all',14);
        keyboard;
    case 'FIG_acute_variability'
        
        subplot(3,3,[7 2])
        rs_imana('acute_variability_all');
        title ('full connectivity pattern')
        
        subplot(3,3,3)
        rs_imana('acute_variability_inter');
        title ('interhemispheric')
        
        subplot(3,3,6)
        rs_imana('acute_variability_intra_les');
        title ('intrahemispheric lesioned')
        
        subplot(3,3,9)
        rs_imana('acute_variability_intra_non');
        title ('intrahemispheric nonlesioned')
        keyboard;
        
    case 'all_distance_patients_longitude_1_4'
        D = load(fullfile(ppDir,'nc_getThemall_jorn.mat'));
        
        %D = getrow(D,D.FM >=40);
        D = getrow(D,D.FM <=40);
        
        %%CAVE
        Dp      = getrow(D,D.control==0);      % change for patients
        
        Dw      = getrow(Dp,Dp.week==1);        % get acute stage
        Dtest   = getrow(Dp,Dp.week==4);        % get week that should be compared
        
        % calculate real Euclidian distance
        D.euc_week_real = pdist2((mean(Dw.M,1)), (mean(Dtest.M,1)),'euclidean');
        
        % permutation
        P = [];
        h = waitbar(0,'Please wait, I am busy Bootstraping the Data');
        NumPer=10000;
        
        % make fake groups by mixing weeks
        Dfake = getrow(Dp,ismember(Dp.week,[1 4]));
        
        SN          = [1:length(Dfake.subj_name)]';
        Dfake.SN    = SN;
        
        for j=1:NumPer
            %snRef       = datasample(SN,11,'Replace',false);
            snRef       = datasample(SN,7,'Replace',false);
            snNonRef    = SN(~ismember(SN,snRef));
            
            m1 = mean(Dfake.M(snRef,:),1);
            m2 = mean(Dfake.M(snNonRef,:),1);
            
            Pi.euc_week_fake              = pdist2(m1,m2,'euclidean');
            
            Pi.j = j;
            
            h=waitbar(j / NumPer);
            P = addstruct(P,Pi);
            
        end
        
        Numb_exceed       = (D.euc_week_real<P.euc_week_fake);
        EUC_per_all       = (sum(Numb_exceed))/10000;
        keyboard
        
        % Histoplot Figures
        % nhist(P.euc_week_fake,'noerror','smooth','color','sequential','samebins');
        % drawline(D.euc_week_real,'dir','vert','linestyle','--','color',[1 0 0])
        %    y = prctile(P.euc_week_fake, [5 95]);
        %    drawline(y,'dir','vert','linestyle',':','color',[0 0 0])
        % set_graphics(gcf,'ylabel',{'frequency'}, 'xlabel',{'Euclidian distance'},'fontsize_all',14);
        
        %2.Figure type
        y = prctile(P.euc_week_fake, [2.5 97.5]);
        x1=0.95;
        x2=1.05;
        y1=y(1);
        y2=y(2);
        y = [y1, y1, y2, y2, y1];
        x = [x1, x2, x2, x1, x1];
        plot(x, y, 'b-', 'LineWidth', 3);
        hold on
        plot(1,D.euc_week_real,'r*')
        xlim([0.9 1.1])
        ylim([0.5 1.5])
    case 'all_distance_patients_longitude_1_12'
        D = load(fullfile(ppDir,'nc_getThemall_jorn.mat'));
        %D = getrow(D,D.FM >=40);
        D = getrow(D,D.FM <=40);
        
        %%CAVE
        Dp      = getrow(D,D.control==0);
        
        Dw      = getrow(Dp,Dp.week==1);
        Dtest   = getrow(Dp,Dp.week==12);
        
        
        D.euc_week_real = pdist2((mean(Dw.M,1)), (mean(Dtest.M,1)),'euclidean');
        
        % permutation
        P = [];
        h = waitbar(0,'Please wait, Bootstraping the Data');
        NumPer=10000;
        
        Dfake = getrow(Dp,ismember(Dp.week,[1 12]));
        
        SN          = [1:length(Dfake.subj_name)]';
        Dfake.SN    = SN;
        
        for j=1:NumPer
            %snRef       = datasample(SN,11,'Replace',false);
            snRef       = datasample(SN,7,'Replace',false);
            snNonRef    = SN(~ismember(SN,snRef));
            
            m1 = mean(Dfake.M(snRef,:),1);
            m2 = mean(Dfake.M(snNonRef,:),1);
            
            Pi.euc_week_fake              = pdist2(m1,m2,'euclidean');
            
            Pi.j = j;
            
            h=waitbar(j / NumPer);
            P = addstruct(P,Pi);
            
        end
        
        Numb_exceed       = (D.euc_week_real<P.euc_week_fake);
        EUC_per_all       = (sum(Numb_exceed))/10000;
        keyboard
        
        % Histoplot Figures
        % nhist(P.euc_week_fake,'noerror','smooth','color','sequential','samebins');
        % drawline(D.euc_week_real,'dir','vert','linestyle','--','color',[1 0 0])
        % y = prctile(P.euc_week_fake, [5 95]);
        % drawline(y,'dir','vert','linestyle',':','color',[0 0 0])
        % set_graphics(gcf,'ylabel',{'frequency'}, 'xlabel',{'Euclidian distance'},'fontsize_all',14);
        
        % 2. Figure
        y = prctile(P.euc_week_fake, [2.5 97.5]);
        x1=1.95;
        x2=2.05;
        y1=y(1);
        y2=y(2);
        y = [y1, y1, y2, y2, y1];
        x = [x1, x2, x2, x1, x1];
        plot(x, y, 'b-', 'LineWidth', 3);
        hold on
        plot(2,D.euc_week_real,'r*')
        xlim([1.9 2.1])
        ylim([0.5 1.5])
    case 'all_distance_patients_longitude_1_24'
        D = load(fullfile(ppDir,'nc_getThemall_jorn.mat'));
        %D = getrow(D,D.FM >=40);
        D = getrow(D,D.FM <=40);
        
        %%CAVE
        Dp      = getrow(D,D.control==0);
        
        Dw      = getrow(Dp,Dp.week==1);
        Dtest   = getrow(Dp,Dp.week==24);
        
        
        D.euc_week_real = pdist2((mean(Dw.M,1)), (mean(Dtest.M,1)),'euclidean');
        
        % permutation
        P = [];
        h = waitbar(0,'Please wait, Bootstraping the Data');
        NumPer=10000;
        
        Dfake = getrow(Dp,ismember(Dp.week,[1 24]));
        
        SN          = [1:length(Dfake.subj_name)]';
        Dfake.SN    = SN;
        
        for j=1:NumPer
            %snRef       = datasample(SN,11,'Replace',false);
            snRef       = datasample(SN,7,'Replace',false);
            snNonRef    = SN(~ismember(SN,snRef));
            
            m1 = mean(Dfake.M(snRef,:),1);
            m2 = mean(Dfake.M(snNonRef,:),1);
            
            Pi.euc_week_fake              = pdist2(m1,m2,'euclidean');
            
            Pi.j = j;
            
            h=waitbar(j / NumPer);
            P = addstruct(P,Pi);
            
        end
        
        Numb_exceed       = (D.euc_week_real<P.euc_week_fake);
        EUC_per_all       = (sum(Numb_exceed))/10000;
        keyboard
        
        %  nhist(P.euc_week_fake,'noerror','smooth','color','sequential','samebins');
        %  drawline(D.euc_week_real,'dir','vert','linestyle','--','color',[1 0 0])
        %  y = prctile(P.euc_week_fake, [5 95]);
        %  drawline(y,'dir','vert','linestyle',':','color',[0 0 0])
        %  set_graphics(gcf,'ylabel',{'frequency'}, 'xlabel',{'Euclidian distance'},'fontsize_all',14);
        
        % 2. Figure
        y = prctile(P.euc_week_fake, [2.5 97.5]);
        x1=2.95;
        x2=3.05;
        y1=y(1);
        y2=y(2);
        y = [y1, y1, y2, y2, y1];
        x = [x1, x2, x2, x1, x1];
        plot(x, y, 'b-', 'LineWidth', 3);
        hold on
        plot(3,D.euc_week_real,'r*')
        xlim([2.9 3.1])
        ylim([0.5 1.5])
    case 'all_distance_patients_longitude_1_52'
        D = load(fullfile(ppDir,'nc_getThemall_jorn.mat'));
        %D = getrow(D,D.FM >=40);
        D = getrow(D,D.FM <=40);
        
        %%CAVE
        Dp      = getrow(D,D.control==0);
        
        Dw      = getrow(Dp,Dp.week==1);
        Dtest   = getrow(Dp,Dp.week==52);
        
        
        D.euc_week_real = pdist2((mean(Dw.M,1)), (mean(Dtest.M,1)),'euclidean');
        
        % permutation
        P = [];
        h = waitbar(0,'Please wait, Bootstraping the Data');
        NumPer=10000;
        
        Dfake = getrow(Dp,ismember(Dp.week,[1 52]));
        
        SN          = [1:length(Dfake.subj_name)]';
        Dfake.SN    = SN;
        
        for j=1:NumPer
            %snRef       = datasample(SN,11,'Replace',false);
            snRef       = datasample(SN,7,'Replace',false);
            snNonRef    = SN(~ismember(SN,snRef));
            
            m1 = mean(Dfake.M(snRef,:),1);
            m2 = mean(Dfake.M(snNonRef,:),1);
            
            Pi.euc_week_fake              = pdist2(m1,m2,'euclidean');
            
            Pi.j = j;
            
            h=waitbar(j / NumPer);
            P = addstruct(P,Pi);
            
        end
        
        Numb_exceed       = (D.euc_week_real<P.euc_week_fake);
        EUC_per_all       = (sum(Numb_exceed))/10000;
        keyboard
        
        % nhist(P.euc_week_fake,'noerror','smooth','color','sequential','samebins');
        % drawline(D.euc_week_real,'dir','vert','linestyle','--','color',[1 0 0])
        % y = prctile(P.euc_week_fake, [5 95]);
        % drawline(y,'dir','vert','linestyle',':','color',[0 0 0])
        % set_graphics(gcf,'ylabel',{'frequency'}, 'xlabel',{'Euclidian distance'},'fontsize_all',14);
        
        % 2. Figure
        y = prctile(P.euc_week_fake, [2.5 97.5]);
        x1=3.95;
        x2=4.05;
        y1=y(1);
        y2=y(2);
        y = [y1, y1, y2, y2, y1];
        x = [x1, x2, x2, x1, x1];
        plot(x, y, 'b-', 'LineWidth', 3);
        hold on
        plot(4,D.euc_week_real,'r*')
        xlim([3.9 4.1])
        ylim([0.5 1.5])
    case 'FIG_all_week_patient'
        subplot(1,4,1)
        rs_imana('all_distance_patients_longitude_1_4');
        title ('All patients week 1 versus 4')
        
        subplot(1,4,2)
        rs_imana('all_distance_patients_longitude_1_12');
        title ('All patients week 1 versus 12')
        
        subplot(1,4,3)
        rs_imana('all_distance_patients_longitude_1_24');
        title ('All patients week 1 versus 24')
        
        subplot(1,4,4)
        rs_imana('all_distance_patients_longitude_1_52');
        title ('All patients week 1 versus 52')
        keyboard;
        
    case 'inter_distance_patients_longitude_1_4'
        D = load(fullfile(ppDir,'nc_interhem.mat'));
        D = getrow(D,D.FM >=40);
        
        Dp      = getrow(D,D.control==0);
        Dw      = getrow(Dp,Dp.week==1);
        Dtest   = getrow(Dp,Dp.week==4);
        
        
        D.euc_week_real = pdist2((mean(Dw.M,1)), (mean(Dtest.M,1)),'euclidean');
        
        % permutation
        P = [];
        h = waitbar(0,'Please wait, Bootstraping the Data');
        NumPer=10000;
        
        Dfake = getrow(Dp,ismember(Dp.week,[1 4]));
        
        SN          = [1:length(Dfake.subj_name)]';
        Dfake.SN    = SN;
        
        for j=1:NumPer
            snRef       = datasample(SN,11,'Replace',false);
            snNonRef    = SN(~ismember(SN,snRef));
            
            m1 = mean(Dfake.M(snRef,:),1);
            m2 = mean(Dfake.M(snNonRef,:),1);
            
            Pi.euc_week_fake              = pdist2(m1,m2,'euclidean');
            
            Pi.j = j;
            
            h=waitbar(j / NumPer);
            P = addstruct(P,Pi);
            
        end
        
        Numb_exceed       = (D.euc_week_real<P.euc_week_fake);
        EUC_per_all       = (sum(Numb_exceed))/10000;
        keyboard
        
        %     nhist(P.euc_week_fake,'noerror','smooth','color','sequential','samebins');
        %     drawline(D.euc_week_real,'dir','vert','linestyle','--','color',[1 0 0])
        %     set_graphics(gcf,'ylabel',{'frequency'}, 'xlabel',{'Euclidian distance'},'fontsize_all',14);
        
        % 2. Figure
        y = prctile(P.euc_week_fake, [2.5 97.5]);
        x1=3.95;
        x2=4.05;
        y1=y(1);
        y2=y(2);
        y = [y1, y1, y2, y2, y1];
        x = [x1, x2, x2, x1, x1];
        plot(x, y, 'b-', 'LineWidth', 3);
        hold on
        plot(4,D.euc_week_real,'r*')
        xlim([3.9 4.1])
        ylim([0.5 1.5])
    case 'inter_distance_patients_longitude_1_12'
        D = load(fullfile(ppDir,'nc_interhem.mat'));
        D = getrow(D,D.FM >=40);
        
        Dp      = getrow(D,D.control==0);
        Dw      = getrow(Dp,Dp.week==1);
        Dtest   = getrow(Dp,Dp.week==12);
        
        
        D.euc_week_real = pdist2((mean(Dw.M,1)), (mean(Dtest.M,1)),'euclidean');
        
        % permutation
        P = [];
        h = waitbar(0,'Please wait, Bootstraping the Data');
        NumPer=10000;
        
        Dfake = getrow(Dp,ismember(Dp.week,[1 12]));
        
        SN          = [1:length(Dfake.subj_name)]';
        Dfake.SN    = SN;
        
        for j=1:NumPer
            snRef       = datasample(SN,11,'Replace',false);
            snNonRef    = SN(~ismember(SN,snRef));
            
            m1 = mean(Dfake.M(snRef,:),1);
            m2 = mean(Dfake.M(snNonRef,:),1);
            
            Pi.euc_week_fake              = pdist2(m1,m2,'euclidean');
            
            Pi.j = j;
            
            h=waitbar(j / NumPer);
            P = addstruct(P,Pi);
            
        end
        
        Numb_exceed       = (D.euc_week_real<P.euc_week_fake);
        EUC_per_all       = (sum(Numb_exceed))/10000;
        keyboard
        
        %     nhist(P.euc_week_fake,'noerror','smooth','color','sequential','samebins');
        %     drawline(D.euc_week_real,'dir','vert','linestyle','--','color',[1 0 0])
        %     set_graphics(gcf,'ylabel',{'frequency'}, 'xlabel',{'Euclidian distance'},'fontsize_all',14);
        
        % 2. Figure
        y = prctile(P.euc_week_fake, [2.5 97.5]);
        x1=3.95;
        x2=4.05;
        y1=y(1);
        y2=y(2);
        y = [y1, y1, y2, y2, y1];
        x = [x1, x2, x2, x1, x1];
        plot(x, y, 'b-', 'LineWidth', 3);
        hold on
        plot(4,D.euc_week_real,'r*')
        xlim([3.9 4.1])
        ylim([0.5 1.5])
    case 'inter_distance_patients_longitude_1_24'
        D = load(fullfile(ppDir,'nc_interhem.mat'));
        D = getrow(D,D.FM >=40);
        
        Dp      = getrow(D,D.control==0);
        Dw      = getrow(Dp,Dp.week==1);
        Dtest   = getrow(Dp,Dp.week==24);
        
        
        D.euc_week_real = pdist2((mean(Dw.M,1)), (mean(Dtest.M,1)),'euclidean');
        
        % permutation
        P = [];
        h = waitbar(0,'Please wait, Bootstraping the Data');
        NumPer=10000;
        
        Dfake = getrow(Dp,ismember(Dp.week,[1 24]));
        
        SN          = [1:length(Dfake.subj_name)]';
        Dfake.SN    = SN;
        
        for j=1:NumPer
            snRef       = datasample(SN,11,'Replace',false);
            snNonRef    = SN(~ismember(SN,snRef));
            
            m1 = mean(Dfake.M(snRef,:),1);
            m2 = mean(Dfake.M(snNonRef,:),1);
            
            Pi.euc_week_fake              = pdist2(m1,m2,'euclidean');
            
            Pi.j = j;
            
            h=waitbar(j / NumPer);
            P = addstruct(P,Pi);
            
        end
        
        Numb_exceed       = (D.euc_week_real<P.euc_week_fake);
        EUC_per_all       = (sum(Numb_exceed))/10000;
        keyboard
        
        %     nhist(P.euc_week_fake,'noerror','smooth','color','sequential','samebins');
        %     drawline(D.euc_week_real,'dir','vert','linestyle','--','color',[1 0 0])
        %     set_graphics(gcf,'ylabel',{'frequency'}, 'xlabel',{'Euclidian distance'},'fontsize_all',14);
        
        % 2. Figure
        y = prctile(P.euc_week_fake, [2.5 97.5]);
        x1=3.95;
        x2=4.05;
        y1=y(1);
        y2=y(2);
        y = [y1, y1, y2, y2, y1];
        x = [x1, x2, x2, x1, x1];
        plot(x, y, 'b-', 'LineWidth', 3);
        hold on
        plot(4,D.euc_week_real,'r*')
        xlim([3.9 4.1])
        ylim([0.5 1.5])
    case 'inter_distance_patients_longitude_1_52'
        D = load(fullfile(ppDir,'nc_interhem.mat'));
        D = getrow(D,D.FM >=40);
        
        Dp      = getrow(D,D.control==0);
        Dw      = getrow(Dp,Dp.week==1);
        Dtest   = getrow(Dp,Dp.week==52);
        
        
        D.euc_week_real = pdist2((mean(Dw.M,1)), (mean(Dtest.M,1)),'euclidean');
        
        % permutation
        P = [];
        h = waitbar(0,'Please wait, Bootstraping the Data');
        NumPer=10000;
        
        Dfake = getrow(Dp,ismember(Dp.week,[1 52]));
        
        SN          = [1:length(Dfake.subj_name)]';
        Dfake.SN    = SN;
        
        for j=1:NumPer
            snRef       = datasample(SN,11,'Replace',false);
            snNonRef    = SN(~ismember(SN,snRef));
            
            m1 = mean(Dfake.M(snRef,:),1);
            m2 = mean(Dfake.M(snNonRef,:),1);
            
            Pi.euc_week_fake              = pdist2(m1,m2,'euclidean');
            
            Pi.j = j;
            
            h=waitbar(j / NumPer);
            P = addstruct(P,Pi);
            
        end
        
        Numb_exceed       = (D.euc_week_real<P.euc_week_fake);
        EUC_per_all       = (sum(Numb_exceed))/10000;
        keyboard
        
        %     nhist(P.euc_week_fake,'noerror','smooth','color','sequential','samebins');
        %     drawline(D.euc_week_real,'dir','vert','linestyle','--','color',[1 0 0])
        %     set_graphics(gcf,'ylabel',{'frequency'}, 'xlabel',{'Euclidian distance'},'fontsize_all',14);
        
        % 2. Figure
        y = prctile(P.euc_week_fake, [2.5 97.5]);
        x1=3.95;
        x2=4.05;
        y1=y(1);
        y2=y(2);
        y = [y1, y1, y2, y2, y1];
        x = [x1, x2, x2, x1, x1];
        plot(x, y, 'b-', 'LineWidth', 3);
        hold on
        plot(4,D.euc_week_real,'r*')
        xlim([3.9 4.1])
        ylim([0.5 1.5])
    case 'FIG_inter_week_patient'
        subplot(1,4,1)
        rs_imana('inter_distance_patients_longitude_1_4');
        title ('inter patients week 1 versus 4')
        
        subplot(1,4,2)
        rs_imana('inter_distance_patients_longitude_1_12');
        title ('inter patients week 1 versus 12')
        
        subplot(1,4,3)
        rs_imana('inter_distance_patients_longitude_1_24');
        title ('inter patients week 1 versus 24')
        
        subplot(1,4,4)
        rs_imana('inter_distance_patients_longitude_1_52');
        title ('inter patients week 1 versus 52')
        keyboard;
        
    case 'intra_non_distance_patients_longitude_1_4'
        D = load(fullfile(ppDir,'nc_intrahem_nonlesioned.mat'));
        D = getrow(D,D.FM >=40);
        
        Dp      = getrow(D,D.control==0);
        Dw      = getrow(Dp,Dp.week==1);
        Dtest   = getrow(Dp,Dp.week==4);
        
        
        D.euc_week_real = pdist2((mean(Dw.M,1)), (mean(Dtest.M,1)),'euclidean');
        
        % permutation
        P = [];
        h = waitbar(0,'Please wait, Bootstraping the Data');
        NumPer=10000;
        
        Dfake = getrow(Dp,ismember(Dp.week,[1 4]));
        
        SN          = [1:length(Dfake.subj_name)]';
        Dfake.SN    = SN;
        
        for j=1:NumPer
            snRef       = datasample(SN,11,'Replace',false);
            snNonRef    = SN(~ismember(SN,snRef));
            
            m1 = mean(Dfake.M(snRef,:),1);
            m2 = mean(Dfake.M(snNonRef,:),1);
            
            Pi.euc_week_fake              = pdist2(m1,m2,'euclidean');
            
            Pi.j = j;
            
            h=waitbar(j / NumPer);
            P = addstruct(P,Pi);
            
        end
        
        Numb_exceed       = (D.euc_week_real<P.euc_week_fake);
        EUC_per_all       = (sum(Numb_exceed))/10000;
        keyboard
        
        %     nhist(P.euc_week_fake,'noerror','smooth','color','sequential','samebins');
        %     drawline(D.euc_week_real,'dir','vert','linestyle','--','color',[1 0 0])
        %     set_graphics(gcf,'ylabel',{'frequency'}, 'xlabel',{'Euclidian distance'},'fontsize_all',14);
        
        % 2. Figure
        y = prctile(P.euc_week_fake, [2.5 97.5]);
        x1=3.95;
        x2=4.05;
        y1=y(1);
        y2=y(2);
        y = [y1, y1, y2, y2, y1];
        x = [x1, x2, x2, x1, x1];
        plot(x, y, 'b-', 'LineWidth', 3);
        hold on
        plot(4,D.euc_week_real,'r*')
        xlim([3.9 4.1])
        ylim([0.5 1.5])
    case 'intra_non_distance_patients_longitude_1_12'
        D = load(fullfile(ppDir,'nc_intrahem_nonlesioned.mat'));
        D = getrow(D,D.FM >=40);
        
        Dp      = getrow(D,D.control==0);
        Dw      = getrow(Dp,Dp.week==1);
        Dtest   = getrow(Dp,Dp.week==12);
        
        
        D.euc_week_real = pdist2((mean(Dw.M,1)), (mean(Dtest.M,1)),'euclidean');
        
        % permutation
        P = [];
        h = waitbar(0,'Please wait, Bootstraping the Data');
        NumPer=10000;
        
        Dfake = getrow(Dp,ismember(Dp.week,[1 12]));
        
        SN          = [1:length(Dfake.subj_name)]';
        Dfake.SN    = SN;
        
        for j=1:NumPer
            snRef       = datasample(SN,11,'Replace',false);
            snNonRef    = SN(~ismember(SN,snRef));
            
            m1 = mean(Dfake.M(snRef,:),1);
            m2 = mean(Dfake.M(snNonRef,:),1);
            
            Pi.euc_week_fake              = pdist2(m1,m2,'euclidean');
            
            Pi.j = j;
            
            h=waitbar(j / NumPer);
            P = addstruct(P,Pi);
            
        end
        
        Numb_exceed       = (D.euc_week_real<P.euc_week_fake);
        EUC_per_all       = (sum(Numb_exceed))/10000;
        keyboard
        
        %     nhist(P.euc_week_fake,'noerror','smooth','color','sequential','samebins');
        %     drawline(D.euc_week_real,'dir','vert','linestyle','--','color',[1 0 0])
        %     set_graphics(gcf,'ylabel',{'frequency'}, 'xlabel',{'Euclidian distance'},'fontsize_all',14);
        
        % 2. Figure
        y = prctile(P.euc_week_fake, [2.5 97.5]);
        x1=3.95;
        x2=4.05;
        y1=y(1);
        y2=y(2);
        y = [y1, y1, y2, y2, y1];
        x = [x1, x2, x2, x1, x1];
        plot(x, y, 'b-', 'LineWidth', 3);
        hold on
        plot(4,D.euc_week_real,'r*')
        xlim([3.9 4.1])
        ylim([0.5 1.5])
    case 'intra_non_distance_patients_longitude_1_24'
        D = load(fullfile(ppDir,'nc_intrahem_nonlesioned.mat'));
        D = getrow(D,D.FM >=40);
        
        Dp      = getrow(D,D.control==0);
        Dw      = getrow(Dp,Dp.week==1);
        Dtest   = getrow(Dp,Dp.week==24);
        
        
        D.euc_week_real = pdist2((mean(Dw.M,1)), (mean(Dtest.M,1)),'euclidean');
        
        % permutation
        P = [];
        h = waitbar(0,'Please wait, Bootstraping the Data');
        NumPer=10000;
        
        Dfake = getrow(Dp,ismember(Dp.week,[1 24]));
        
        SN          = [1:length(Dfake.subj_name)]';
        Dfake.SN    = SN;
        
        for j=1:NumPer
            snRef       = datasample(SN,11,'Replace',false);
            snNonRef    = SN(~ismember(SN,snRef));
            
            m1 = mean(Dfake.M(snRef,:),1);
            m2 = mean(Dfake.M(snNonRef,:),1);
            
            Pi.euc_week_fake              = pdist2(m1,m2,'euclidean');
            
            Pi.j = j;
            
            h=waitbar(j / NumPer);
            P = addstruct(P,Pi);
            
        end
        
        Numb_exceed       = (D.euc_week_real<P.euc_week_fake);
        EUC_per_all       = (sum(Numb_exceed))/10000;
        keyboard
        
        %     nhist(P.euc_week_fake,'noerror','smooth','color','sequential','samebins');
        %     drawline(D.euc_week_real,'dir','vert','linestyle','--','color',[1 0 0])
        %     set_graphics(gcf,'ylabel',{'frequency'}, 'xlabel',{'Euclidian distance'},'fontsize_all',14);
        
        % 2. Figure
        y = prctile(P.euc_week_fake, [2.5 97.5]);
        x1=3.95;
        x2=4.05;
        y1=y(1);
        y2=y(2);
        y = [y1, y1, y2, y2, y1];
        x = [x1, x2, x2, x1, x1];
        plot(x, y, 'b-', 'LineWidth', 3);
        hold on
        plot(4,D.euc_week_real,'r*')
        xlim([3.9 4.1])
        ylim([0.5 1.5])
    case 'intra_non_distance_patients_longitude_1_52'
        D = load(fullfile(ppDir,'nc_intrahem_nonlesioned.mat'));
        D = getrow(D,D.FM >=40);
        
        Dp      = getrow(D,D.control==0);
        Dw      = getrow(Dp,Dp.week==1);
        Dtest   = getrow(Dp,Dp.week==52);
        
        
        D.euc_week_real = pdist2((mean(Dw.M,1)), (mean(Dtest.M,1)),'euclidean');
        
        % permutation
        P = [];
        h = waitbar(0,'Please wait, Bootstraping the Data');
        NumPer=10000;
        
        Dfake = getrow(Dp,ismember(Dp.week,[1 52]));
        
        SN          = [1:length(Dfake.subj_name)]';
        Dfake.SN    = SN;
        
        for j=1:NumPer
            snRef       = datasample(SN,11,'Replace',false);
            snNonRef    = SN(~ismember(SN,snRef));
            
            m1 = mean(Dfake.M(snRef,:),1);
            m2 = mean(Dfake.M(snNonRef,:),1);
            
            Pi.euc_week_fake              = pdist2(m1,m2,'euclidean');
            
            Pi.j = j;
            
            h=waitbar(j / NumPer);
            P = addstruct(P,Pi);
            
        end
        
        Numb_exceed       = (D.euc_week_real<P.euc_week_fake);
        EUC_per_all       = (sum(Numb_exceed))/10000;
        keyboard
        
        %     nhist(P.euc_week_fake,'noerror','smooth','color','sequential','samebins');
        %     drawline(D.euc_week_real,'dir','vert','linestyle','--','color',[1 0 0])
        %     set_graphics(gcf,'ylabel',{'frequency'}, 'xlabel',{'Euclidian distance'},'fontsize_all',14);
        
        % 2. Figure
        y = prctile(P.euc_week_fake, [2.5 97.5]);
        x1=3.95;
        x2=4.05;
        y1=y(1);
        y2=y(2);
        y = [y1, y1, y2, y2, y1];
        x = [x1, x2, x2, x1, x1];
        plot(x, y, 'b-', 'LineWidth', 3);
        hold on
        plot(4,D.euc_week_real,'r*')
        xlim([3.9 4.1])
        ylim([0.5 1.5])
    case 'FIG_intra_non_week_patient'
        subplot(1,4,1)
        rs_imana('intra_non_distance_patients_longitude_1_4');
        title ('intra non patients week 1 versus 4')
        
        subplot(1,4,2)
        rs_imana('intra_non_distance_patients_longitude_1_12');
        title ('intra non patients week 1 versus 12')
        
        subplot(1,4,3)
        rs_imana('intra_non_distance_patients_longitude_1_24');
        title ('intra non patients week 1 versus 24')
        
        subplot(1,4,4)
        rs_imana('intra_non_distance_patients_longitude_1_52');
        title ('intra non patients week 1 versus 52')
        keyboard;
        
    case 'intra_les_distance_patients_longitude_1_4'
        D = load(fullfile(ppDir,'nc_intrahem_lesioned.mat'));
        D = getrow(D,D.FM >=40);
        
        Dp      = getrow(D,D.control==0);
        Dw      = getrow(Dp,Dp.week==1);
        Dtest   = getrow(Dp,Dp.week==4);
        
        
        D.euc_week_real = pdist2((mean(Dw.M,1)), (mean(Dtest.M,1)),'euclidean');
        
        % permutation
        P = [];
        h = waitbar(0,'Please wait, Bootstraping the Data');
        NumPer=10000;
        
        Dfake = getrow(Dp,ismember(Dp.week,[1 4]));
        
        SN          = [1:length(Dfake.subj_name)]';
        Dfake.SN    = SN;
        
        for j=1:NumPer
            snRef       = datasample(SN,11,'Replace',false);
            snNonRef    = SN(~ismember(SN,snRef));
            
            m1 = mean(Dfake.M(snRef,:),1);
            m2 = mean(Dfake.M(snNonRef,:),1);
            
            Pi.euc_week_fake              = pdist2(m1,m2,'euclidean');
            
            Pi.j = j;
            
            h=waitbar(j / NumPer);
            P = addstruct(P,Pi);
            
        end
        
        Numb_exceed       = (D.euc_week_real<P.euc_week_fake);
        EUC_per_all       = (sum(Numb_exceed))/10000;
        keyboard
        
        %     nhist(P.euc_week_fake,'noerror','smooth','color','sequential','samebins');
        %     drawline(D.euc_week_real,'dir','vert','linestyle','--','color',[1 0 0])
        %     set_graphics(gcf,'ylabel',{'frequency'}, 'xlabel',{'Euclidian distance'},'fontsize_all',14);
        
        % 2. Figure
        y = prctile(P.euc_week_fake, [2.5 97.5]);
        x1=3.95;
        x2=4.05;
        y1=y(1);
        y2=y(2);
        y = [y1, y1, y2, y2, y1];
        x = [x1, x2, x2, x1, x1];
        plot(x, y, 'b-', 'LineWidth', 3);
        hold on
        plot(4,D.euc_week_real,'r*')
        xlim([3.9 4.1])
        ylim([0.5 1.5])
    case 'intra_les_distance_patients_longitude_1_12'
        D = load(fullfile(ppDir,'nc_intrahem_lesioned.mat'));
        D = getrow(D,D.FM >=40);
        
        Dp      = getrow(D,D.control==0);
        Dw      = getrow(Dp,Dp.week==1);
        Dtest   = getrow(Dp,Dp.week==12);
        
        
        D.euc_week_real = pdist2((mean(Dw.M,1)), (mean(Dtest.M,1)),'euclidean');
        
        % permutation
        P = [];
        h = waitbar(0,'Please wait, Bootstraping the Data');
        NumPer=10000;
        
        Dfake = getrow(Dp,ismember(Dp.week,[1 12]));
        
        SN          = [1:length(Dfake.subj_name)]';
        Dfake.SN    = SN;
        
        for j=1:NumPer
            snRef       = datasample(SN,11,'Replace',false);
            snNonRef    = SN(~ismember(SN,snRef));
            
            m1 = mean(Dfake.M(snRef,:),1);
            m2 = mean(Dfake.M(snNonRef,:),1);
            
            Pi.euc_week_fake              = pdist2(m1,m2,'euclidean');
            
            Pi.j = j;
            
            h=waitbar(j / NumPer);
            P = addstruct(P,Pi);
            
        end
        
        Numb_exceed       = (D.euc_week_real<P.euc_week_fake);
        EUC_per_all       = (sum(Numb_exceed))/10000;
        keyboard
        
        %     nhist(P.euc_week_fake,'noerror','smooth','color','sequential','samebins');
        %     drawline(D.euc_week_real,'dir','vert','linestyle','--','color',[1 0 0])
        %     set_graphics(gcf,'ylabel',{'frequency'}, 'xlabel',{'Euclidian distance'},'fontsize_all',14);
        
        % 2. Figure
        y = prctile(P.euc_week_fake, [2.5 97.5]);
        x1=3.95;
        x2=4.05;
        y1=y(1);
        y2=y(2);
        y = [y1, y1, y2, y2, y1];
        x = [x1, x2, x2, x1, x1];
        plot(x, y, 'b-', 'LineWidth', 3);
        hold on
        plot(4,D.euc_week_real,'r*')
        xlim([3.9 4.1])
        ylim([0.5 1.5])
    case 'intra_les_distance_patients_longitude_1_24'
        D = load(fullfile(ppDir,'nc_intrahem_lesioned.mat'));
        D = getrow(D,D.FM >=40);
        
        Dp      = getrow(D,D.control==0);
        Dw      = getrow(Dp,Dp.week==1);
        Dtest   = getrow(Dp,Dp.week==24);
        
        
        D.euc_week_real = pdist2((mean(Dw.M,1)), (mean(Dtest.M,1)),'euclidean');
        
        % permutation
        P = [];
        h = waitbar(0,'Please wait, Bootstraping the Data');
        NumPer=10000;
        
        Dfake = getrow(Dp,ismember(Dp.week,[1 24]));
        
        SN          = [1:length(Dfake.subj_name)]';
        Dfake.SN    = SN;
        
        for j=1:NumPer
            snRef       = datasample(SN,11,'Replace',false);
            snNonRef    = SN(~ismember(SN,snRef));
            
            m1 = mean(Dfake.M(snRef,:),1);
            m2 = mean(Dfake.M(snNonRef,:),1);
            
            Pi.euc_week_fake              = pdist2(m1,m2,'euclidean');
            
            Pi.j = j;
            
            h=waitbar(j / NumPer);
            P = addstruct(P,Pi);
            
        end
        
        Numb_exceed       = (D.euc_week_real<P.euc_week_fake);
        EUC_per_all       = (sum(Numb_exceed))/10000;
        keyboard
        
        %     nhist(P.euc_week_fake,'noerror','smooth','color','sequential','samebins');
        %     drawline(D.euc_week_real,'dir','vert','linestyle','--','color',[1 0 0])
        %     set_graphics(gcf,'ylabel',{'frequency'}, 'xlabel',{'Euclidian distance'},'fontsize_all',14);
        
        % 2. Figure
        y = prctile(P.euc_week_fake, [2.5 97.5]);
        x1=3.95;
        x2=4.05;
        y1=y(1);
        y2=y(2);
        y = [y1, y1, y2, y2, y1];
        x = [x1, x2, x2, x1, x1];
        plot(x, y, 'b-', 'LineWidth', 3);
        hold on
        plot(4,D.euc_week_real,'r*')
        xlim([3.9 4.1])
        ylim([0.5 1.5])
    case 'intra_les_distance_patients_longitude_1_52'
        D = load(fullfile(ppDir,'nc_intrahem_lesioned.mat'));
        D = getrow(D,D.FM >=40);
        
        Dp      = getrow(D,D.control==0);
        Dw      = getrow(Dp,Dp.week==1);
        Dtest   = getrow(Dp,Dp.week==52);
        
        
        D.euc_week_real = pdist2((mean(Dw.M,1)), (mean(Dtest.M,1)),'euclidean');
        
        % permutation
        P = [];
        h = waitbar(0,'Please wait, Bootstraping the Data');
        NumPer=10000;
        
        Dfake = getrow(Dp,ismember(Dp.week,[1 52]));
        
        SN          = [1:length(Dfake.subj_name)]';
        Dfake.SN    = SN;
        
        for j=1:NumPer
            snRef       = datasample(SN,11,'Replace',false);
            snNonRef    = SN(~ismember(SN,snRef));
            
            m1 = mean(Dfake.M(snRef,:),1);
            m2 = mean(Dfake.M(snNonRef,:),1);
            
            Pi.euc_week_fake              = pdist2(m1,m2,'euclidean');
            
            Pi.j = j;
            
            h=waitbar(j / NumPer);
            P = addstruct(P,Pi);
            
        end
        
        Numb_exceed       = (D.euc_week_real<P.euc_week_fake);
        EUC_per_all       = (sum(Numb_exceed))/10000;
        keyboard
        
        %     nhist(P.euc_week_fake,'noerror','smooth','color','sequential','samebins');
        %     drawline(D.euc_week_real,'dir','vert','linestyle','--','color',[1 0 0])
        %     set_graphics(gcf,'ylabel',{'frequency'}, 'xlabel',{'Euclidian distance'},'fontsize_all',14);
        
        % 2. Figure
        y = prctile(P.euc_week_fake, [2.5 97.5]);
        x1=3.95;
        x2=4.05;
        y1=y(1);
        y2=y(2);
        y = [y1, y1, y2, y2, y1];
        x = [x1, x2, x2, x1, x1];
        plot(x, y, 'b-', 'LineWidth', 3);
        hold on
        plot(4,D.euc_week_real,'r*')
        xlim([3.9 4.1])
        ylim([0.5 1.5])
    case 'FIG_intra_les_week_patient'
        subplot(1,4,1)
        rs_imana('intra_les_distance_patients_longitude_1_4');
        title ('intra les patients week 1 versus 4')
        
        subplot(1,4,2)
        rs_imana('intra_les_distance_patients_longitude_1_12');
        title ('intra les patients week 1 versus 12')
        
        subplot(1,4,3)
        rs_imana('intra_les_distance_patients_longitude_1_24');
        title ('intra les patients week 1 versus 24')
        
        subplot(1,4,4)
        rs_imana('intra_les_distance_patients_longitude_1_52');
        title ('intra les patients week 1 versus 52')
        keyboard;
        
    case 'longitudinal_variability_all'                     %bc it is a subset with only people that have all time points different analysis
        D = load(fullfile(ppDir,'nc_getThemall_jorn.mat'));
        
        Ref     = 1;
        Comp    = [4, 12, 24, 52];
        
        % get guys that have measurements for all timepoints
        Dp      = getrow(D,D.control==0);
        [x,sn]  = pivottable(Dp.subj_name,Dp.week,Dp.week,'length');
        sn      = sn(~sum(isnan(x),2));
        Dp      = getrow(Dp,ismember(Dp.subj_name,sn));
        
        % loop over all comparison time points to calculate diff with ref
        S = [];
        for w=Comp
            % loop over each subjects
            for s=sn'
                m1 = getrow(Dp,strcmp(Dp.subj_name,s) & Dp.week==Ref);
                m2 = getrow(Dp,strcmp(Dp.subj_name,s) & Dp.week==w);
                
                d = pdist2(m1.M, m2.M,'euclidean');
                
                Si.sn        = s;
                Si.week      = w;
                Si.euc_dist  = d;
                S = addstruct(S,Si);
            end;
        end;
        
        % stats for anova
        x = pivottable(S.sn,S.week,S.euc_dist,'mean');
        anova1(x);
        
        W1_4 = mean(x(:,1));
        std1_4 = std(x(:,1));
        W1_12 = mean(x(:,2));
        std1_12 = std(x(:,2));
        W1_24 = mean(x(:,3));
        std1_24 = std(x(:,3));
        W1_52 = mean(x(:,4));
        std1_52 = std(x(:,4));
        keyboard;
    case 'longitudinal_variability_inter'
        D = load(fullfile(ppDir,'nc_interhem.mat'));
        
        Ref     = 1;
        Comp    = [4, 12, 24, 52];
        
        % get guys that have measurements for all timepoints
        Dp      = getrow(D,D.control==0);
        [x,sn]  = pivottable(Dp.subj_name,Dp.week,Dp.week,'length');
        sn      = sn(~sum(isnan(x),2));
        Dp      = getrow(Dp,ismember(Dp.subj_name,sn));
        
        % loop over all comparison time points to calculate diff with ref
        S = [];
        for w=Comp
            % loop over each subjects
            for s=sn'
                m1 = getrow(Dp,strcmp(Dp.subj_name,s) & Dp.week==Ref);
                m2 = getrow(Dp,strcmp(Dp.subj_name,s) & Dp.week==w);
                
                d = pdist2(m1.M, m2.M,'euclidean');
                
                Si.sn        = s;
                Si.week      = w;
                Si.euc_dist  = d;
                S = addstruct(S,Si);
            end;
        end;
        
        % stats for anova
        x = pivottable(S.sn,S.week,S.euc_dist,'mean');
        anova1(x);
        
        W1_4 = mean(x(:,1));
        std1_4 = std(x(:,1));
        W1_12 = mean(x(:,2));
        std1_12 = std(x(:,2));
        W1_24 = mean(x(:,3));
        std1_24 = std(x(:,3));
        W1_52 = mean(x(:,4));
        std1_52 = std(x(:,4));
        keyboard;
    case 'longitudinal_variability_intra_non'
        D = load(fullfile(ppDir,'nc_intrahem_nonlesioned.mat'));
        
        Ref     = 1;
        Comp    = [4, 12, 24, 52];
        
        % get guys that have measurements for all timepoints
        Dp      = getrow(D,D.control==0);
        [x,sn]  = pivottable(Dp.subj_name,Dp.week,Dp.week,'length');
        sn      = sn(~sum(isnan(x),2));
        Dp      = getrow(Dp,ismember(Dp.subj_name,sn));
        
        % loop over all comparison time points to calculate diff with ref
        S = [];
        for w=Comp
            % loop over each subjects
            for s=sn'
                m1 = getrow(Dp,strcmp(Dp.subj_name,s) & Dp.week==Ref);
                m2 = getrow(Dp,strcmp(Dp.subj_name,s) & Dp.week==w);
                
                d = pdist2(m1.M, m2.M,'euclidean');
                
                Si.sn        = s;
                Si.week      = w;
                Si.euc_dist  = d;
                S = addstruct(S,Si);
            end;
        end;
        
        % stats for anova
        x = pivottable(S.sn,S.week,S.euc_dist,'mean');
        anova1(x);
        
        W1_4 = mean(x(:,1));
        std1_4 = std(x(:,1));
        W1_12 = mean(x(:,2));
        std1_12 = std(x(:,2));
        W1_24 = mean(x(:,3));
        std1_24 = std(x(:,3));
        W1_52 = mean(x(:,4));
        std1_52 = std(x(:,4));
        keyboard;
    case 'longitudinal_variability_intra_les'
        D = load(fullfile(ppDir,'nc_intrahem_lesioned.mat'));
        
        Ref     = 1;
        Comp    = [4, 12, 24, 52];
        
        % get guys that have measurements for all timepoints
        Dp      = getrow(D,D.control==0);
        [x,sn]  = pivottable(Dp.subj_name,Dp.week,Dp.week,'length');
        sn      = sn(~sum(isnan(x),2));
        Dp      = getrow(Dp,ismember(Dp.subj_name,sn));
        
        % loop over all comparison time points to calculate diff with ref
        S = [];
        for w=Comp
            % loop over each subjects
            for s=sn'
                m1 = getrow(Dp,strcmp(Dp.subj_name,s) & Dp.week==Ref);
                m2 = getrow(Dp,strcmp(Dp.subj_name,s) & Dp.week==w);
                
                d = pdist2(m1.M, m2.M,'euclidean');
                
                Si.sn        = s;
                Si.week      = w;
                Si.euc_dist  = d;
                S = addstruct(S,Si);
            end;
        end;
        
        % stats for anova
        x = pivottable(S.sn,S.week,S.euc_dist,'mean');
        anova1(x);
        
        W1_4 = mean(x(:,1));
        std1_4 = std(x(:,1));
        W1_12 = mean(x(:,2));
        std1_12 = std(x(:,2));
        W1_24 = mean(x(:,3));
        std1_24 = std(x(:,3));
        W1_52 = mean(x(:,4));
        std1_52 = std(x(:,4));
        keyboard;
        
    case 'plot_networkCorr_inter'
        % interhemispheric: plots heatmaps, vectors, correlation
        D   = load(fullfile(ppDir,'nc_interhem.mat'));
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
        traceplot(1:25,fisherinv((Dc.M)),'errorfcn','stderr');
        hold on
        traceplot(1:25,fisherinv((DpW1.M)),'errorfcn','stderr','linecolor',[1 0 0],'patchcolor',[1 0 0]);
        ylim([0 1]);        xlim([0 25]);
        subplot(232);
        traceplot(1:25,fisherinv((Dc.M)),'errorfcn','stderr');
        hold on
        traceplot(1:25,fisherinv((DpW4.M)),'errorfcn','stderr','linecolor',[1 0 0],'patchcolor',[1 0 0]);
        ylim([0 1]);        xlim([0 25]);
        title('Week4');
        ylabel('correlation index')
        xlabel('roi pairs')
        subplot(233);
        traceplot(1:25,fisherinv((Dc.M)),'errorfcn','stderr');
        hold on
        traceplot(1:25,fisherinv((DpW12.M)),'errorfcn','stderr','linecolor',[1 0 0],'patchcolor',[1 0 0]);
        title('Week12');
        ylabel('correlation index')
        xlabel('roi pairs')
        ylim([0 1]);        xlim([0 25]);
        subplot(234);
        traceplot(1:25,fisherinv((Dc.M)),'errorfcn','stderr');
        hold on
        traceplot(1:25,fisherinv((DpW24.M)),'errorfcn','stderr','linecolor',[1 0 0],'patchcolor',[1 0 0]);
        title('Week24');
        ylabel('correlation index')
        xlabel('roi pairs')
        ylim([0 1]);        xlim([0 25]);
        subplot(235);
        traceplot(1:25,fisherinv((Dc.M)),'errorfcn','stderr');
        hold on
        traceplot(1:25,fisherinv((DpW52.M)),'errorfcn','stderr','linecolor',[1 0 0],'patchcolor',[1 0 0]);
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
    case 'plot_all_conn'
        D = load(fullfile(ppDir,'nc_getThemall_jorn.mat'));
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
        colormapeditor
        
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
    case 'plot_networkCorr_lesioned'    % plot controls and patients (lesioned hemisphere only)
        D   = load(fullfile(ppDir,'rs_networkCorr_lesionedIntra.mat'));
        
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
    case 'STATS_behavioural_timecourse'
        D = rs_imana('PP_studySubset');
        D.mvc = nanmean(D.mvc,2);
        
        D =     tapply(D,{'subj_name','week','control'},{'FM','nanmean'},...
            {'ARAT','nanmean'},{'mvc','nanmean'});
        dsave(fullfile(statsDir,'behavioural_timecourse.dat'),D);
        
        x = pivottable(D.subj_name,D.week,D.FM,'nanmean','subset',D.control==1);
        y = pivottable(D.subj_name,D.week,D.FM,'nanmean','subset',D.control==0);
        ttest(x(:,1),y(:,1), 2,'independent');
        
        x = pivottable(D.subj_name,D.week,D.ARAT,'nanmean','subset',D.control==1);
        y = pivottable(D.subj_name,D.week,D.ARAT,'nanmean','subset',D.control==0);
        ttest(x(:,1),y(:,1), 2,'independent');
        
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
    case 'ROI_volumes'
        load (fullfile(roiDir,'UZ_2365_regions_all.mat'));
        for i=[1 2 3 4 5]
            fprintf('Regio %s has %d voxels\n',R{i}.name,size(R{i}.data,1));
        end;
    otherwise
        disp('no such case');
end;