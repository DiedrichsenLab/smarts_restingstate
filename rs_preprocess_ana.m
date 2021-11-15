function varargout = rs_preprocess_ana(what,varargin)
% what - the type of data processing to be performed

% Directories

rootDir         = '/Volumes/porsche';
rootDir         = '/Volumes/diedrichsen_data$';

%rsDir           = fullfile(rootDir,'/meret/MATLAB/data_rsfMRI/resting_state_for_Jorn');

rsDir           = fullfile(rootDir,'/data/smarts/fmri/restingstate_imaging/');
rssDir          = fullfile(rootDir,'/data/smarts/fmri/restingstate_imaging/preprocess/files');
anaDir          = fullfile(rootDir,'anatomicals');


regName = {'S1','M1','PMd','PMv','SMA','V1','IPS','OPJ'};
regSide = [1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2]';
regType = [1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8]';
weeks   = [1 4 12 24 52];
Weeks   = {'W1','W4','W12','W24','W52'};

numVol      = 200;
nDiscard    = 10;
hem     = {'lh','rh'};
hemName = {'LeftHem','RightHem'};


switch(what)
    case 'LIST'
        fname = fullfile(rsDir,'patient_list.txt');

        if nargout==0
            x = readtable(fname,'Delimiter','\t','ReadRowNames',1);
            disp(x);
        else
            D           = dload(fname);
            D.subj_name = strcat(D.Centre,'_',num2str(D.ID));
            varargout = {D};
        end;
    case 'check_data_exists'
        D   = rs_preprocess_ana('LIST');

        T = [];
        for i=1:length(D.SN)
            Di = getrow(D,i);
            Ti = Di;
            
            %f = {'meanWIPBoldRestSENSE_001.img','preproc_Bold_Rest_'};
            f = {'meanBoldRestSENSE_001.img','preproc_Bold_Rest_'};
            
            ppDir     = fullfile(rssDir);
            rawDir    = fullfile(rsDir,Di.subj_name{1},Di.Week{1});
            rawFile = fullfile(rawDir,f{1});
            ppFile  = fullfile(ppDir,[f{2} Di.subj_name{1} '_' Di.Week{1} '.nii']);
    
            Ti.existRawImg  = logical(exist(rawFile));
            Ti.existPPImg   = logical(exist(ppFile));
            
            if nargout == 0
                fprintf('%d.\t%s\t\traw:%d\tpp:%d\n',Di.SN,[Di.subj_name{1} ' ' Di.Week{1}],...
                            Ti.existRawImg,Ti.existPPImg);
            end;
            
            T       = addstruct(T,Ti);
        end;        
        
        if nargout > 0
            varargout = {T};
        end;
    case 'make_mean_nii'
        % 0. check which subjects have the raw and preprocessed files
        D = rs_preprocess_ana('check_data_exists');
        D = getrow(D,D.existPPImg & D.existRawImg);
            
        % 1. fix problems in variables/file names etc
        %   - fix problem in one control where week is not specified
        %   (default to W0)
        idx             = cellfun(@length,D.RefT1) == 1;
        D.RefT1(idx)    = {'W0'};

        % 2. set up automatic co-registration
        %   - create copy of mean epi
        %   - bias correct mean epi
        %   - align mean epi to anatomical
        for s = 1:length(D.RefT1)
            Di = getrow(D,s);
            
            f = fullfile(rsDir,Di.subj_name{1},Di.Week{1},'meanBoldRestSENSE_001.hdr');
            v = spm_vol(f);
            x = spm_read_vols(v);
            
            v.fname = fullfile(rsDir,Di.subj_name{1},Di.Week{1},'rs_meanepi.nii');
            spm_write_vol(v,x);
            
            fprintf('%d\n',s);
        end;
    case 'bias_correct'
        % 0. check which subjects have the raw and preprocessed files
        D = rs_preprocess_ana('check_data_exists');
        D = getrow(D,D.existPPImg & D.existRawImg);
            
        % 1. fix problems in variables/file names etc
        %   - fix problem in one control where week is not specified
        %   (default to W0)
        idx             = cellfun(@length,D.RefT1) == 1;
        D.RefT1(idx)    = {'W0'};

        % 2. set up automatic co-registration
        %   - create copy of mean epi
        %   - bias correct mean epi
        %   - align mean epi to anatomical
        for s = 1:length(D.RefT1)
            Di = getrow(D,s);
            
            f = fullfile(rsDir,Di.subj_name{1},Di.Week{1},'rs_meanepi.nii');
            spmj_bias_correct({f});
            
            fprintf('%d\n',s);
        end;
    case 'coreg'
        % 0. check which subjects have the raw and preprocessed files
        D = rs_preprocess_ana('check_data_exists');
        D = getrow(D,D.existPPImg & D.existRawImg);
            
        % 1. fix problems in variables/file names etc
        %   - fix problem in one control where week is not specified
        %   (default to W0)
        idx             = cellfun(@length,D.RefT1) == 1;
        D.RefT1(idx)    = {'W0'};

        % 2. set up automatic co-registration
        %   - create copy of mean epi
        %   - bias correct mean epi
        %   - align mean epi to anatomical
        for s = 1:length(D.RefT1)
            Di = getrow(D,s);
            
            f   = fullfile(rsDir,Di.subj_name{1},Di.Week{1},'brs_meanepi.nii');
            ana = fullfile(anaDir,Di.subj_name{1},Di.RefT1{1},[Di.subj_name{1},'_',Di.RefT1{1},'_T1.nii']);
            
            J.ref = {ana};
            J.source = {f};
            J.other = {''};
            J.eoptions.cost_fun = 'nmi';
            J.eoptions.sep = [4 2];
            J.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
            J.eoptions.fwhm = [7 7];
            matlabbatch{1}.spm.spatial.coreg.estimate=J;
            spm_jobman('run',matlabbatch);
            
            fprintf('%d\n',s);
        end;
    %% change for different preprocessing pipelines
    case 'make_same_align'

        % 0. check which subjects have the raw and preprocessed files
        D = rs_preprocess_ana('check_data_exists');
        D = getrow(D,D.existPPImg & D.existRawImg);
            
        % 1. fix problems in variables/file names etc
        %   - fix problem in one control where week is not specified
        %   (default to W0)
        idx             = cellfun(@length,D.RefT1) == 1;
        D.RefT1(idx)    = {'W0'};

        % 2. set up automatic co-registration
        %   - create copy of mean epi
        %   - bias correct mean epi
        %   - align mean epi to anatomical
        for s = 1:length(D.RefT1)
            Di = getrow(D,s);
            
            %ftype = {'pre_Bold_Rest','preproc_Bold_Rest'};
            ftype = {'Denoised_Bold_Rest','preproc_Bold_Rest','pre_Bold_Rest'};         %after reviewers requested pipeline

            % reference mean epi which contains affine transform matrix to copy
            P{1} = fullfile(rsDir,Di.subj_name{1},Di.Week{1},'brs_meanepi.nii');
            
            for i = 1:2
                f = fullfile(rsDir,Di.subj_name{1},Di.Week{1},[ftype{i},'_',Di.subj_name{1},'_',Di.Week{1},'.nii']);
                
                if exist(f)
                    Q = {};
                    for j=1:numVol
                        Q{end+1}    = sprintf('%s,%d',f,j);
                    end;
                    
                    % align individual resting state images
                    spmj_makesamealign_nifti(char(P),char(Q));
                    
                    fprintf('%d. %s, %d\n',s,Di.subj_name{1},i);
                end;
            end;
        end;
    case 'check_same_align'
        % 0. check which subjects have the raw and preprocessed files
        D = rs_preprocess_ana('check_data_exists');
        D = getrow(D,D.existPPImg & D.existRawImg);
            
        % 1. fix problems in variables/file names etc
        %   - fix problem in one control where week is not specified
        %   (default to W0)
        idx             = cellfun(@length,D.RefT1) == 1;
        D.RefT1(idx)    = {'W0'};

        % 2. set up automatic co-registration
        %   - create copy of mean epi
        %   - bias correct mean epi
        %   - align mean epi to anatomical
        for s = 1:length(D.RefT1)
            Di = getrow(D,s);
            
            %ftype = {'pre_Bold_Rest','preproc_Bold_Rest'};
            ftype = {'Denoised_Bold_Rest','preproc_Bold_Rest','pre_Bold_Rest'};

            % reference mean epi which contains affine transform matrix to copy
            P{1} = fullfile(rsDir,Di.subj_name{1},Di.Week{1},'brs_meanepi.nii');
            
            for i = 1:2
                f = fullfile(rsDir,Di.subj_name{1},Di.Week{1},[ftype{i},'_',Di.subj_name{1},'_',Di.Week{1},'.nii']);
                
                if exist(f)
                    Q = {};
                    for j=1:numVol
                        Q{end+1}    = sprintf('%s,%d',f,j);
                    end;
                    
                    % align individual resting state images
                    spmj_checksamealign(char(P),char(Q));
                    
                    fprintf('%d. %s, %d\n',s,Di.subj_name{1},i);
                end;
            end;
        end;
                
end;
