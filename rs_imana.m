function varargout=rs_imana(what,varargin)
% MB - analysis script for the SMARTS resting state data

% rootDir             = '/Volumes/External/data/smarts';
% rootDir             = '/Users/naveed/Documents/data/smarts';
% rootDir             = '/Users/meret/Dropbox/N&M/rsfMRI_results';
% rootDir             = '/Users/naveed/Dropbox/Personal/N&M/rsfMRI_results';
rootDir             = '/Volumes/Naveed/data/smarts';
behDir              = [rootDir '/bedside/analysis'];
rsDir               = [rootDir '/fmri/restingstate_imaging'];
roiDir              = [rootDir '/fmri/RegionOfInterest'];
ppDir               = [rsDir '/preprocess/'];

% color_scheme
colours             = num2cell(colormap(parula(20)),2)';
sty_2grp_ac         = colours([8,12]);   % e.g. group 1 vs group 2
sty_multiple        = colours([10;14;16]);


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
    case 'process_all'
        pre = {'preproc_Bold_Rest'};
        lesionType = {'none','cortical','subcortical','mixed'}; %1,2,3,4
        %roi  = {[1 9],[2 10],[3 11],[4 12],[5 13],[6 14],[7 15],[8 16],... %{'S1','M1','PMd','PMv','SMA','V1','IPS','OPJ'};
        %[2 1],[2 3],[2 4],[2 5],[2 6],[10 11],[10 12],[10 13],[10 14],[10 15]};                                                                 %1-8 is left, 9-16 is right
%         roi = [combntns(1:17,2)]'; % all possible pairs not including OJP
        roi  = combnk(1:16,2); % all pairs
        
        cd(ppDir);
        d = dir('*.nii');
        
        D   = dload(fullfile(ppDir,'patient_list.txt'));
        B   = load(fullfile(behDir,'ens_alldat.mat'));
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
                % Si.M1_involv  = percentage of M1 lesion involvement
                % Si.CST_involv = percentage of CST lesion involvement
                if isempty(find(strcmp(D.LesionLocation(i),lesionType)))
                    Si.lesionType = nan;
                end;
                Si.lesionSide   = D.lesionside(i);
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
                    Si.r        = rp;
                    S           = addstruct(S,Si);
                end;
            end;
        end;
        varargout = {S};
        save(fullfile(ppDir,'rs_preprocess.mat'),'-struct','S');        
    case 'cf_connectivity'
    D   = load(fullfile(rootDir,'rs_preprocess.mat'));  
    %only take subcorticals for now
    %or define threshold
    D   = getrow(D,(D.lesionType~=2)); %no cortical
    D   = getrow(D,(D.lesionType~=4)); %no mixed
    
%    Needs to be updated when all roi pairs are included
%     M1  = getrow(D,D.rType==2); 
%     Pmd = getrow(D,D.rType==3);
%     Pmv = getrow(D,D.rType==4);
%     V1  = getrow(D,D.rType==6);
%     
%     M_SMA   = getrow(D,D.rType==11);
%     M_Pmd   = getrow(D,D.rType==9);
%     M_Pmv   = getrow(D,D.rType==10);
%     M_V1    = getrow(D,D.rType==6);
    
    s   = style_sheet(sty_2grp_ac,'leg',{'patients','control'},'leglocation','northoutside');
     
    %interhemispheric connectivity
    h = make_figure;
    subplot(2,2,1)
    title('M1-M1 subcortical')
    lineplot(M1.week,M1.r,'plotfcn','nanmean','errorfcn','stderr','split',M1.control,'CAT',s(1).CAT,s(1).PLOT{:});
    ylabel('correlation index')
    xlabel('weeks')
    subplot(2,2,2)
    title('Pmd-Pmd subcortical')
    lineplot(Pmd.week,Pmd.r,'plotfcn','nanmean','errorfcn','stderr','split',M1.control,'CAT',s(1).CAT,s(1).PLOT{:});
    ylabel('correlation index')
    xlabel('weeks')
    subplot(2,2,3)
    title('Pmv-Pmv subcortical')
    lineplot(Pmv.week,Pmv.r,'plotfcn','nanmean','errorfcn','stderr','split',M1.control,'CAT',s(1).CAT,s(1).PLOT{:});
    ylabel('correlation index')
    xlabel('weeks')
    subplot(2,2,4)
    title('V1-V1 subcortical')
    lineplot(V1.week,V1.r,'plotfcn','nanmean','errorfcn','stderr','split',M1.control,'CAT',s(1).CAT,s(1).PLOT{:});
    ylabel('correlation index')
    xlabel('weeks') 
    
    % intracortical connectivity
    % how to choose only lesioned side?
%     h = make_figure;
%     subplot(4,4,1)
%     title('M1-SMA subcortical')
%     lineplot(M1.week,M1.r,'plotfcn','nanmean','errorfcn','stderr','split',M1.control,'CAT',s(1).CAT,s(1).PLOT{:});
%     ylabel('correlation index')
%     xlabel('weeks')
%     subplot(4,4,2)
%     title('M1-Pmd subcortical')
%     lineplot(Pmd.week,Pmd.r,'plotfcn','nanmean','errorfcn','stderr','split',M1.control,'CAT',s(1).CAT,s(1).PLOT{:});
%     ylabel('correlation index')
%     xlabel('weeks')
%     subplot(4,4,3)
%     title('M1-Pmv subcortical')
%     lineplot(Pmv.week,Pmv.r,'plotfcn','nanmean','errorfcn','stderr','split',M1.control,'CAT',s(1).CAT,s(1).PLOT{:});
%     ylabel('correlation index')
%     xlabel('weeks')
%     subplot(2,2,4)
%     title('M1-V1 subcortical')
%     lineplot(V1.week,V1.r,'plotfcn','nanmean','errorfcn','stderr','split',M1.control,'CAT',s(1).CAT,s(1).PLOT{:});
%     ylabel('correlation index')
%     xlabel('weeks')


    case'cf_conn_lesion'
    D   = load(fullfile(rootDir,'rs_preprocess.mat'));  
    D   = load(fullfile(rootDir,'rs_preprocess.mat'));  
    %only take subcorticals for now
    %or define threshold
    D   = getrow(D,(D.lesionType~=2)); %no cortical
    D   = getrow(D,(D.lesionType~=4)); %no mixed
      
%    Needs to be updated when all roi pairs are included
%     M1  = getrow(D,D.rType==2); 

    s   = style_sheet(sty_2grp_ac,'leg',{'patients','control'},'leglocation','northoutside');
     
    %interhemispheric connectivity
    h = make_figure;
    subplot(2,2,1)
    scatterplot(M1.CST_involve,M1.r,'regression','linear','printcorr','split',M1.control,'markercolor',sty_2grp_ac);
    
    case'cf_conn_behav'
    D   = load(fullfile(rootDir,'rs_preprocess.mat'));  
    
    % only take subcorticals for now
    D   = getrow(D,(D.lesionType~=2)); %no cortical
    D   = getrow(D,(D.lesionType~=4)); %no mixed
    
    case 'make_networkCorr'         % correlate motor networks for each subj/week
        D   = load(fullfile(rootDir,'rs_preprocess.mat')); 
        S   = [];
        for sn=unique(D.subj_name)'
            Ds = getrow(D,strcmp(D.subj_name,sn));
            for w=unique(Ds.week)'
                Dw  = getrow(Ds,Ds.week==w);
                fprintf('%s - W%d\n',Dw.subj_name{1},Dw.week(1));
                
                M   = zeros(16,16);     % matrix of corr. for all regions
                for i=1:size(Dw.rpair,1)
                    M(Dw.rpair(i,1),Dw.rpair(i,2)) = Dw.r(i);
                end;
                M = M+M';       % symmetricizing matrix
                
                % add data to new structure
                Si      = getrow(Dw,1);
                Si      = rmfield(Si,{'rType','rpair','r'});
                Si.M    = nonsym_squareform(M);
                S       = addstruct(S,Si);
            end;
        end;
        save(fullfile(rootDir,'rs_networkCorr.mat'),'-struct','S');                
    case 'plot_networkCorr'     % plot correlation across controls and patients
        D   = load(fullfile(rootDir,'rs_networkCorr.mat'));                
        Dc  = getrow(D,D.control==1);
        Dp  = getrow(D,D.control~=1);
        
        % get avg matrices for each group
        Mc  = nanmean(Dc.M,1);
        Mp  = nanmean(Dp.M,1);
        Mc  = nonsym_squareform(Mc);
        Mp  = nonsym_squareform(Mp);
        
        scale = [0 0.8];
        h=make_figure;
        subplot(121);
        colormap(hot); 
        imagesc(Mc,scale); colorbar; axis square;
        title('Controls');
        subplot(122);
        colormap(hot); 
        imagesc(Mp,scale); colorbar; axis square;
        title('Patients');
        set_graphics(h,'xtick',1:16,'ytick',1:16)
        keyboard;
    otherwise
        disp('no such case');
end;