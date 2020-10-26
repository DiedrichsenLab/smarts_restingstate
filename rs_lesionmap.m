function varargout = rs_lesionmap (what,varargin)


rootDir     = '/Volumes/Naveed/data/smarts/fmri';
rsDir       = fullfile(rootDir,'/restingstate_imaging/preprocess');
anaDir      = fullfile(rootDir,'anatomicals');
lesionDir   = fullfile(rootDir,'MNI_space');


subj_id  = {'CUP_1002','CU_2310','CU_2663','CU_2925','JHP_1001','JHP_1002',...
            'JHP_1004','JHU_2395','JHU_2531','JHU_2789','JHU_3176','UZP_1001',...
            'UZP_1002','UZP_1004','UZP_1005','UZP_1006','UZP_1007','UZP_1008',...
            'UZ_2365','UZ_2450','UZ_2565','UZ_2652','UZ_2654','UZ_3239','UZ_3240',...
            'UZ_3241','UZ_3243','UZ_3246','UZ_3247','UZ_3248'};
        
switch(what)
    case 'get_subj'
        D = load(fullfile(rsDir,'rs_preprocess.mat'));
        D = getrow(D,ismember(D.subj_name,subj_id));
        
        varargout = {D};
        
    case 'average_lesion'
        S = rs_lesionmap('get_subj');
        S = getrow(S,S.control==0); 
        S = tapply(S,{'subj_name','week'},{'lesionType','mean'},{'lesionSide','mean'});

        % get subjs
        sn = unique(S.subj_name)';
        
        % go into lesion dir (lesions have been mapped to eve/mni space so that all 
        % in the same space)
        cd(fullfile(anaDir,'MNI_space')); 
        
        L = []; % matrix of lesion image file names
       
        % loop over all subjects
        for s=1:length(sn)
            indx = find(strcmp(S.subj_name,sn{s}));     % Find all sessions from this subject
            si   = indx(1);                             % Take first session - assumes that it always orderd in text file
            L    = spm_vol(['w' S.subj_name{si} '_lesion.nii']); 
            
            X(:,:,:,s)=spm_read_vols(L); 
            
            % flip left lesions to right
            if S.lesionSide(si)==1 % left
                disp(['flipping lesion from left to right for ' S.subj_name{si}]);
                X(:,:,:,s) = flipdim(X(:,:,:,s),1);
            end;
            
            % write flipped image
            [~,fname] = fileparts(L.fname);
            L.fname   = [fname '_flipped.nii'];
            spm_write_vol(L,X(:,:,:,s)>0);
        end; 
        Avrg=mean(X>0,4); 
        
        L.fname = fullfile('restingstate_lesion_overlap_mni.nii');
        L.dt= [16 0]; % Data type set to float 
        spm_write_vol(L,Avrg);
end;