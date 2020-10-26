function varargout = premovement_analysis(what,varargin)
% function varargout = premovement_analysis(what,varargin)
% analysis procedures for analyzing SMARTS premovement TMS data
% Jing Xu 02/25/2016; changed Meret 03/20/2016

exp_code = 'IHI'; % 'SICI'

baseDir = '/Users/meret/Documents/MATLAB/preprocess-selected';
dataDir = fullfile(baseDir,'preprocess_debug',exp_code);
analysisDir = fullfile(baseDir,'Analysis_debug',exp_code);
if strcmp(exp_code,'IHI')
    subject_file = 'IHI_list.txt';
elseif strcmp(exp_code,'SICI')
    subject_file = 'SICI_list.txt';
end    

S = dload(fullfile(dataDir,subject_file));
nSubj = length(S.SN);
for s=1:nSubj
    S.code1{s,1} = [S.Centre{s} num2str(S.ID(s)) '_' S.Week{s}]; % Unique session name, use for paradigm data
    S.code2{s,1} = [S.Centre{s} num2str(S.ID(s)) '.' S.Muscle{s} S.LesionSide{s} '.' S.Week{s}]; % Unique session name, used for Signal data
    S.name{s,1}  = [S.Centre{s} num2str(S.ID(s))];
end;
% for weeks
weeks = [1;4;12;24;52];
color_w={[0 0 0],[0.5 0 0],[1 0 0],[1 0.5 0],[1 1 0]};
Weeks={'W1','W4','W12','W24','W52'};

switch (what)
    case 'subject_list'              % Prints subject / session list
        if nargin < 2
            sn = [1:nSubj];
        else
            sn = varargin{1};
        end;
        
        fprintf('\nSN\tCentre\tID\tWeek\tweek\tblocks\tLesionSide\tlesionside\tMuscle\tisPatient');
        fprintf('\n--\t------\t----\t----\t----\t------\t-----------\t----------\t-----\t----\n');
        for s = 1:nSubj
            % printing subject lists
            fprintf([num2str(S.SN(s)) '\t' S.Centre{s} '\t' num2str(S.ID(s)) '\t' S.Week{s} ... 
                '\t' num2str(S.week(s)) '\t' num2str(S.blocks(s)) '\t' S.LesionSide{s}...
                '\t\t' num2str(S.lesionside(s)) '\t\t' S.Muscle{s} '\t' num2str(S.isPatient(s)) '\n']);
        end;
        fprintf('\n');
        varargout = {S};        
    case 'compare_markings'          % Helper case used in QA: comparing two people's preprocess markings
        dir1 = '/Volumes/dayu/McDonnell/TMS/data/preprocess/IHI/Meret';
        dir2 = '/Volumes/dayu/McDonnell/TMS/data/preprocess/IHI';
        if nargin < 2
            sn = [1:nSubj];
        else
            sn = varargin{1};
        end;
        for s=sn
            for b=1:S.blocks(s)
                fname = [exp_code(1:3) 'P' num2str(b) '.' S.code2{s} '.MEP.mat'];
                disp(['Comparing ' fname]);
                clear D1 D2;
                if exist(fullfile(dir1,S.name{s},S.Week{s},fname),'file')
                    D1 = load(fullfile(dir1,S.name{s},S.Week{s},fname));
                else
                    disp(['File does not exist in ' dir1]);
                end
                if exist(fullfile(dir2,S.name{s},S.Week{s},fname),'file')
                    D2 = load(fullfile(dir2,S.name{s},S.Week{s},fname));
                else
                    disp(['File does not exist in ' dir2]);
                end
                if exist('D1','var')&&exist('D2','var')
                    for ch=1:2
                        Dd(b).channel(ch).number = D1.Summary.channel(ch).number;
                        Dd(b).channel(ch).start  = D1.Summary.channel(ch).start;
                        Dd(b).channel(ch).state  = D1.Summary.channel(ch).state;
                        Dd(b).channel(ch).Dpeak       = D1.Summary.channel(ch).peak - D2.Summary.channel(ch).peak;
                        Dd(b).channel(ch).DpeakTime   = D1.Summary.channel(ch).peakTime - D2.Summary.channel(ch).peakTime;
                        Dd(b).channel(ch).Dtrough     = D1.Summary.channel(ch).trough - D2.Summary.channel(ch).trough;
                        Dd(b).channel(ch).DtroughTime = D1.Summary.channel(ch).troughTime - D2.Summary.channel(ch).troughTime;
                        Dd(b).channel(ch).DRT         = D1.Summary.channel(ch).RT - D2.Summary.channel(ch).RT;
                        Dd(b).channel(ch).Drms        = D1.Summary.channel(ch).rms - D2.Summary.channel(ch).rms;
                        Dd(b).channel(ch).DisGood     = D1.Summary.channel(ch).isGood - D2.Summary.channel(ch).isGood;
                        Dd(b).channel(ch).Damplitude  = D1.Summary.channel(ch).amplitude - D2.Summary.channel(ch).amplitude;
                        Dd(b).channel(ch).Dpre_act_trial_avg = D1.Summary.channel(ch).pre_act_trial_avg - D2.Summary.channel(ch).pre_act_trial_avg;
                        Dd(b).channel(ch).Dpre_act_trial_std = D1.Summary.channel(ch).pre_act_trial_std - D2.Summary.channel(ch).pre_act_trial_std;

                        Dd(b).channel(ch).Dpeak(isnan(D1.Summary.channel(ch).peak))          = D2.Summary.channel(ch).peak(isnan(D1.Summary.channel(ch).peak));
                        Dd(b).channel(ch).Dpeak(isnan(D2.Summary.channel(ch).peak))          = D1.Summary.channel(ch).peak(isnan(D2.Summary.channel(ch).peak));
                        Dd(b).channel(ch).DpeakTime(isnan(D1.Summary.channel(ch).peakTime))  = D2.Summary.channel(ch).peakTime(isnan(D1.Summary.channel(ch).peakTime));
                        Dd(b).channel(ch).DpeakTime(isnan(D2.Summary.channel(ch).peakTime))  = D1.Summary.channel(ch).peakTime(isnan(D2.Summary.channel(ch).peakTime));
                        Dd(b).channel(ch).Dtrough(isnan(D1.Summary.channel(ch).trough))      = D2.Summary.channel(ch).trough(isnan(D1.Summary.channel(ch).trough));
                        Dd(b).channel(ch).Dtrough(isnan(D2.Summary.channel(ch).trough))      = D1.Summary.channel(ch).trough(isnan(D2.Summary.channel(ch).trough));
                        Dd(b).channel(ch).DpeakTime(isnan(D1.Summary.channel(ch).troughTime))= D2.Summary.channel(ch).troughTime(isnan(D1.Summary.channel(ch).troughTime));
                        Dd(b).channel(ch).DpeakTime(isnan(D2.Summary.channel(ch).troughTime))= D1.Summary.channel(ch).troughTime(isnan(D2.Summary.channel(ch).troughTime));
                        Dd(b).channel(ch).DRT(isnan(D1.Summary.channel(ch).RT))              = D2.Summary.channel(ch).RT(isnan(D1.Summary.channel(ch).RT));
                        Dd(b).channel(ch).DRT(isnan(D2.Summary.channel(ch).RT))              = D1.Summary.channel(ch).RT(isnan(D2.Summary.channel(ch).RT));
                        Dd(b).channel(ch).Damplitude(isnan(D1.Summary.channel(ch).amplitude))= D2.Summary.channel(ch).amplitude(isnan(D1.Summary.channel(ch).amplitude));
                        Dd(b).channel(ch).Damplitude(isnan(D2.Summary.channel(ch).amplitude))= D1.Summary.channel(ch).amplitude(isnan(D2.Summary.channel(ch).amplitude));

                        Dd(b).channel(ch).mismatch_amplitude = find(Dd(b).channel(ch).Damplitude~=0&~isnan(Dd(b).channel(ch).Damplitude)...
                                                                    &(D1.Summary.channel(ch).isGood==1 | D2.Summary.channel(ch).isGood==1));
                        Dd(b).channel(ch).mismatch_isGood    = find(Dd(b).channel(ch).DisGood~=0);
                        Dd(b).channel(ch).mismatch_RT        = find(abs(Dd(b).channel(ch).DRT)>10&~isnan(Dd(b).channel(ch).DRT)...
                                                                    &(D1.Summary.channel(ch).isGood==1 | D2.Summary.channel(ch).isGood==1));
                        disp(['B' num2str(b) ' ch' num2str(ch) ': amplitude: ' num2str(Dd(b).channel(ch).mismatch_amplitude')]);
                        disp(['B' num2str(b) ' ch' num2str(ch) ': isGood   : ' num2str(Dd(b).channel(ch).mismatch_isGood')]);
                        disp(['B' num2str(b) ' ch' num2str(ch) ': RT       : ' num2str(Dd(b).channel(ch).mismatch_RT')]);
                    end
                end
            end
        end
        
        varargout = {Dd}; 
    case 'preprocess'                % Gets each subject's data summary
        if nargin < 2
            sn = [1:nSubj];
        else
            sn = varargin{1};
        end;
        for i = sn
            disp(['Processing session no: ' num2str(i)]);
            sessDir = fullfile(dataDir,S.name{i},S.Week{i});
            task = load(fullfile(sessDir,[exp_code '_' S.code1{i} '_paradigm.mat']));
            cd(sessDir);
            
            % check go-cue & tms pulses are presented at the right timing
            task.tms.good_stim = ((task.tms.tms_time-task.tms.go_onset)*100 -...
                                  (task.tms.scheduled*task.tms.pre_gort))<10;
            mep         = []; 
            isGood      = [];
            RT          = []; 
            BN          = [];
            mep_time    = [];
            rms         = [];
            pre_act     = [];
            pre_act_std = [];
            for b = 1:S.blocks(i)
                load(fullfile(sessDir,[exp_code(1:3) 'P' num2str(b) '.' S.code2{i} '.MEP.mat']));
                mep         = [mep;Summary.channel(1).amplitude];
                % using manual markings to check good MEPs for now
                % need to add the criterion of (70+/-30)ms before RT for TS MEPs
                % for secondary look
                isGood      = [isGood;(Summary.channel(1).isGood & Summary.channel(2).isGood)];
                RT          = [RT; Summary.channel(1).RT];
                BN          = [BN; repmat(b,length(Summary.channel(1).RT),1)];
                mep_time    = [mep_time; min(Summary.channel(1).peakTime,Summary.channel(1).troughTime)];
                rms         = [rms; Summary.channel(1).rms];
                pre_act     = [pre_act; Summary.channel(1).pre_act_trial_avg];
                pre_act_std = [pre_act_std; Summary.channel(1).pre_act_trial_avg];
            end
            task.tms.TN        = task.tms.TN(1:length(mep));
            task.tms.states    = task.tms.states(1:length(mep));
            task.tms.scheduled = task.tms.scheduled(1:length(mep));
            task.tms.good_stim = task.tms.good_stim(1:length(mep));
            
            D.BN             = BN;
            D.TN             = task.tms.TN;
            D.control        = repmat(~S.isPatient(i),length(D.BN),1);
            D.lesionside     = repmat(S.lesionside(i),length(D.BN),1);
            D.subj_name      = repmat({S.name{i}},length(D.BN),1);
            D.good_stim      = task.tms.good_stim;
            D.isGood         = (task.tms.good_stim & isGood);
            D.mep            = mep;
            D.pre_gort       = repmat(task.tms.pre_gort,length(D.BN),1); % Go RT tested before TMS session
            D.ts_time        = (task.tms.scheduled).*task.tms.pre_gort+10; % test pulse relative to Go cue onset
            D.mep_time       = D.ts_time+(mep_time-150); % relative to Go cue onset
            D.RT             = D.ts_time+(RT-150);
            D.rms            = rms;
            D.pre_act        = pre_act;
            D.pre_act_std    = pre_act_std;
            D.state          = task.tms.states;
            D.tms_schedule   = task.tms.scheduled;  % TMS timing as scheduled
            D.tms_schedule_r = D.ts_time./D.RT;     % TMS timing in the real trial           
            D.CS_intensity   = repmat(S.CS_intensity(i),length(D.BN),1);
            D.TS_intensity   = repmat(S.TS_intensity(i),length(D.BN),1);
            D.aCST           = repmat(S.aCST(i),length(D.BN),1);
            D.rCST           = repmat(S.rCST(i),length(D.BN),1);
            outfilename = fullfile(analysisDir,[exp_code '_' S.code1{i} '_ana.mat']);
            save(outfilename,'-struct','D');
        end
    case 'preprocess_combined'       % Gets each subject's data summary
        if nargin < 2
            sn = [1:nSubj];
        else
            sn = varargin{1};
        end;
        for i = sn
            disp(['Processing session no: ' num2str(i)]);
            sessDir = fullfile(dataDir,S.name{i},S.Week{i});
            task = load(fullfile(sessDir,[exp_code '_' S.code1{i} '_paradigm.mat']));
            cd(sessDir);
            
            % Code for reading resting TMS measures
            baselineFile = dir('IHIR*MEP.mat');
            if ~isempty(baselineFile)
                bFile = load(baselineFile(1).name);
            else
                disp('No resting baseline file for subject');
            end;
           bFile.Summary.pre_act = bFile.Summary.pre_act_trial_avg;
           bFile.Summary.pre_act_std = bFile.Summary.pre_act_trial_std;
            C=[];
            D=[];
%             C=rmfield(bFile.Summary,{'peak','peakTime','trough','troughTime','start','tag','sweeps','pre_act_trial_avg','pre_act_trial_std'});
            C.state = bFile.Summary.state;                        
            C.rms = bFile.Summary.rms;
            C.mep = bFile.Summary.amplitude;
            C.isGood = bFile.Summary.isGood;
            C.BN = 7*ones(length(C.state),1);
            C.subj_name = repmat(S.name(i),length(C.state),1);
            C.control = repmat(~S.isPatient(i),length(C.BN),1);
            C.CS_intensity_rest   = repmat(S.CS_intensity_rest(i),length(C.BN),1);
            C.TS_intensity_rest   = repmat(S.TS_intensity_rest(i),length(C.BN),1);
            % check go-cue & tms pulses are presented at the right timing
            task.tms.good_stim = ((task.tms.tms_time-task.tms.go_onset)*100 -...
                                  (task.tms.scheduled*task.tms.pre_gort))<10;
            mep         = []; 
            isGood      = [];
            RT          = []; 
            BN          = [];
            mep_time    = [];
            rms         = [];
            pre_act     = [];
            pre_act_std = [];
            for b = 1:S.blocks(i)
                load(fullfile(sessDir,[exp_code(1:3) 'P' num2str(b) '.' S.code2{i} '.MEP.mat']));
                mep         = [mep;Summary.channel(1).amplitude];
                % using manual markings to check good MEPs for now
                % need to add the criterion of (70+/-30)ms before RT for TS MEPs
                % for secondary look
                isGood      = [isGood;(Summary.channel(1).isGood & Summary.channel(2).isGood)];
                RT          = [RT; Summary.channel(1).RT];
                BN          = [BN; repmat(b,length(Summary.channel(1).RT),1)];
                mep_time    = [mep_time; min(Summary.channel(1).peakTime,Summary.channel(1).troughTime)];
                rms         = [rms; Summary.channel(1).rms];
                pre_act     = [pre_act; Summary.channel(1).pre_act_trial_avg];
                pre_act_std = [pre_act_std; Summary.channel(1).pre_act_trial_avg];
            end
            task.tms.TN        = task.tms.TN(1:length(mep));
            task.tms.states    = task.tms.states(1:length(mep));
            task.tms.scheduled = task.tms.scheduled(1:length(mep));
            task.tms.good_stim = task.tms.good_stim(1:length(mep));
            
            D.BN             = BN;
            D.TN             = task.tms.TN;
            D.control        = repmat(~S.isPatient(i),length(D.BN),1);
            D.lesionside     = repmat(S.lesionside(i),length(D.BN),1);
            D.subj_name      = repmat({S.name{i}},length(D.BN),1);
            D.good_stim      = task.tms.good_stim;
            D.isGood         = (task.tms.good_stim & isGood);
            D.mep            = mep;
            D.pre_gort       = repmat(task.tms.pre_gort,length(D.BN),1); % Go RT tested before TMS session
            D.ts_time        = (task.tms.scheduled).*task.tms.pre_gort+10; % test pulse relative to Go cue onset
            D.mep_time       = D.ts_time+(mep_time-150); % relative to Go cue onset
            D.RT             = D.ts_time+(RT-150);
            D.rms            = rms;
            D.pre_act        = pre_act;
            D.pre_act_std    = pre_act_std;
            D.state          = task.tms.states;
            D.tms_schedule   = task.tms.scheduled;  % TMS timing as scheduled
            D.tms_schedule_r = D.ts_time./D.RT;     % TMS timing in the real trial
            D.CS_intensity   = repmat(S.CS_intensity(i),length(D.BN),1);
            D.TS_intensity   = repmat(S.TS_intensity(i),length(D.BN),1);
            D.aCST           = repmat(S.aCST(i),length(D.BN),1);
            D.rCST           = repmat(S.rCST(i),length(D.BN),1);
            D = addstruct(D,C,'row','force');
            outfilename = fullfile(analysisDir,[exp_code '_' S.code1{i} '_ana.mat']);
            save(outfilename,'-struct','D');
        end                      
    case 'make_alldat'               % Makes the alldata structure: Do this before calling other tabs
        allDat = [];
        
        if (~isempty(varargin))
            sn=varargin{1};
        else
            sn=[1:max(S.SN)]; % take all preprocessed data
        end;
        
        for i = sn
            disp(['Processing no: ' num2str(i)]);
            D = load(fullfile(analysisDir,[exp_code '_' S.code1{i} '_ana.mat']));
            
            vec=ones(size(D.BN,1),1);
            D.SN   = S.SN(i)*vec;   % assign session number
            D.ID   = S.ID(i)*vec;   % Assign subject number
            D.week = S.week(i)*vec; % Week of testing
            
            allDat = addstruct(allDat,D);
        end
        % Assigning a unique ID for each subject
        [~,~,allDat.SubjN] = unique(allDat.subj_name,'stable');
        allDatP = getrow(allDat,allDat.control==0);
        allDatC = getrow(allDat,allDat.control==1);
        save(fullfile(analysisDir,'alldat.mat'),'-struct','allDat');
        dsave(fullfile(analysisDir,'alldat.dat'),allDat);
        save(fullfile(analysisDir,'alldat_patient.mat'),'-struct','allDatP');
        save(fullfile(analysisDir,'alldat_control.mat'),'-struct','allDatC');

    case 'make_summary_tms_scheduled'% Averages individual subject's data by 4 scheduled TMS timings
        D = load('alldat.mat');
        D = getrow(D,D.isGood==1);
        T = tapply(D,{'ID','subj_name','SubjN','week','control','lesionside','tms_schedule','state'},...
                   {'pre_gort'},...
                   {'mep'},...
                   {'RT'},...
                   {D.mep,'length(x)','name','mep_count'});
        varargout = {T};
        save(fullfile(analysisDir,'summary_tms_scheduled.mat'),'-struct','T');
        dsave(fullfile(analysisDir,'summary_tms_scheduled.dat'),T);
    case 'make_summary_real_timing'  % Averages individual subject's data by real TMS timings binned to severl categories
        D = load('alldat.mat');
        D = getrow(D,D.good_stim&~(~D.isGood&((D.RT>0&D.pre_act>0.1)|D.RT>1000)));
        % binning the real timing: 0-0.3, 0.3-0.5, 0.5-0.6, 0.6-0.95, into RT
        D.tms_bin(D.tms_schedule_r>0&D.tms_schedule_r<0.3,1) = 1;
        D.tms_bin(D.tms_schedule_r>=0.3&D.tms_schedule_r<0.5,1) = 2;
        D.tms_bin(D.tms_schedule_r>=0.5&D.tms_schedule_r<0.6,1) = 3;
        D.tms_bin(D.tms_schedule_r>=0.6&D.tms_schedule_r<0.95,1)= 4;
        D.tms_bin(D.tms_schedule_r<0|D.tms_schedule_r>=0.95,1)  = 5;
        D = getrow(D,D.tms_bin~=0);
        T = tapply(D,{'ID','subj_name','SubjN','week','control','lesionside','tms_bin','state'},...
                   {'pre_gort'},...
                   {'mep'},...
                   {'RT'},...
                   {D.mep,'length(x)','name','mep_count'});
        varargout = {T};
        save(fullfile(analysisDir,'summary_tms_realtiming.mat'),'-struct','T');
        dsave(fullfile(analysisDir,'summary_tms_realtiming.dat'),T);
        
    case 'plot_tms_timing'           % Plots mean of a varable (e.g. mep, RT) for single vs. conditioned pulses at each TMS timing
        var = varargin{1};
        count   = 9;
        control = 0;
        timing_type = 'scheduled';  % 'scheduled' vs. 'real'
        vararginoptions(varargin(2:end),{'count','control','timing_type'});
        if strcmp(timing_type,'scheduled')
            D = load(fullfile(analysisDir,'summary_tms_scheduled.mat'));
            timing_var = 'tms_schedule';
        elseif strcmp(timing_type,'real')
            D = load(fullfile(analysisDir,'summary_tms_realtiming.mat'));
            timing_var = 'tms_bin';
        end
        
        DD = getrow(D,D.mep_count>=count & D.control==control);
        for w=1:length(weeks)
            [f1,r1,c1]=pivottable(DD.SubjN,DD.(timing_var),DD.(var),'nanmean','subset',DD.week==weeks(w)&DD.state==1);
            [f2,r2,c2]=pivottable(DD.SubjN,DD.(timing_var),DD.(var),'nanmean','subset',DD.week==weeks(w)&DD.state==2);
            [ff1,rr1,cc1]=pivottable(DD.SubjN,DD.control,DD.pre_gort,'nanmean','subset',DD.week==weeks(w));
            % check number of subjects in each TMS timing
            nlabel = [];
            for t=1:length(c1)
                n1(t) = length(find(~isnan(f1(:,t))));
                n2(t) = length(find(~isnan(f2(:,t))));
                nlabel = [nlabel,{['N=' num2str(n1(t))]},{['N=' num2str(n2(t))]}];
            end
            figure(w);
            legs = {'un-conditioned','conditioned'};
            xcoord = barplot(DD.(timing_var),DD.(var),'split',DD.state,'subset',DD.week==weeks(w),...
                'leg',legs);
            if strcmp(var,'RT');
                drawline(nanmean(ff1),'dir','horz','color','r');  % pre_gort
            end
            if strcmp(timing_type,'scheduled')
                if strcmp(var,'RT')
                    if control ==0
                        ylim([0 600]);
                    else
                        ylim([0 400]);
                    end
                elseif strcmp(var,'mep')
                    if control ==0
                        ylim([0 2.5]);
                    else
                        ylim([0 3.5]);
                    end
                end
            elseif strcmp(timing_type,'real')
                if strcmp(var,'RT')
                    if control ==0
                        ylim([0 500]);
                    else
                        ylim([0 400]);
                    end
                elseif strcmp(var,'mep')
                    if control ==0
                        ylim([0 4.5]);
                    else
                        ylim([0 5.5]);
                    end
                end
            end
            yl = ylim;
            text(xcoord-0.4,repmat(yl(2)/8,size(xcoord)),nlabel);
            xlabel('TMS-timing','fontsize',12);
            ylabel(var,'fontsize',12);
            if control==0
                title(['Patient ' Weeks{w}],'fontsize',14);
            else
                title(['Control ' Weeks{w}],'fontsize',14);
            end
        end
        
        keyboard;
    case 'plot_IHI_curve'            % plots IHI curve with diffferent TMS timings at specified time point
        count = 9;
        timing_type = 'scheduled';  % 'scheduled' vs. 'real'
        vararginoptions(varargin,{'count','timing_type'});
        DF = [];
        if strcmp(timing_type,'scheduled')
            D = load(fullfile(analysisDir,'summary_tms_scheduled.mat'));
            timing_var = 'tms_schedule';
        elseif strcmp(timing_type,'real')
            D = load(fullfile(analysisDir,'summary_tms_realtiming.mat'));
            timing_var = 'tms_bin';
            xtlabel={'0-0.3','0.3-0.5','0.5-0.6','0.6-95','in RT'};
        end
        IDs = unique(D.subj_name);
        tms = unique(D.(timing_var));
        for i=1:length(IDs)
            for w=1:length(weeks)
                for t=1:length(tms)
                    SESS = getrow(D,strcmp(D.subj_name,IDs{i}) &...
                                  D.(timing_var)==tms(t) & D.week==weeks(w));
                    if ~isempty(SESS.ID)
                        % un-conditioned test pulse
                        UT_mep       = SESS.mep(SESS.state==1);
                        UT_mep_count = SESS.mep_count(SESS.state==1);
                        % conditioned test pulse
                        CT_mep       = SESS.mep(SESS.state==2);
                        CT_mep_count = SESS.mep_count(SESS.state==2);
                        
                        if UT_mep_count>=count & CT_mep_count>=count
                            TT.perc        = (CT_mep-UT_mep)/UT_mep;
                            TT.ratio       = CT_mep/UT_mep;
                            TT.ID          = SESS.ID(1);
                            TT.subj_name   = {SESS.subj_name{1}};
                            TT.SubjN       = SESS.SubjN(1);
                            TT.week        = SESS.week(1);
                            TT.control     = SESS.control(1);
                            TT.lesionside  = SESS.lesionside(1);
                            TT.(timing_var)= SESS.(timing_var)(1);
                            TT.RT          = nanmean(SESS.RT);
                            DF = addstruct(DF,TT);
                        end
                    end
                end
            end
        end
        keyboard;
        DF.splitvar = DF.(timing_var);
        DP = getrow(DF,DF.control==0);
        DC = getrow(DF,DF.control==1);
        keyboard;
        % mixed-effect model estimate
        RDP = premovement_analysis('mixed_model',DP,DP.perc);
        
        % average over all controls across weeks
        RDC = premovement_analysis('mixed_model',DC,DC.perc);
        
        figure;
        lineplot(RDP.splitF,RDP.y,'split',RDP.week,...
            'errorval',sqrt(RDP.yvar),'errorcap',0,'errorwidth',1.5,...
            'leg',Weeks,'markersize',8,'markercolor',color_w,'markerfill',color_w,...
            'linecolor',color_w,'errorcolor',color_w,'linewidth',1.5);
        drawline(0,'dir','horz');
        set(gca,'ylim',[-0.6 0.4]);
        ylabel('(CT-UT)/UT','fontsize',12);
        xlabel('TMS timing','fontsize',12);
        title('Patient: mixed effect model estimate','fontsize',14);
        if strcmp(timing_type,'real')
            set(gca,'XTickLabel',xtlabel);
        end

        figure;
        lineplot(RDC.splitF,RDC.y,'split',RDC.week,...
            'errorval',sqrt(RDC.yvar),'errorcap',0,'errorwidth',1.5,...
            'leg',Weeks,'markersize',8,'markercolor',color_w,'markerfill',color_w,...
            'linecolor',color_w,'errorcolor',color_w,'linewidth',1.5);
        drawline(0,'dir','horz');
        set(gca,'ylim',[-0.6 0.4]);
        ylabel('(CT-UT)/UT','fontsize',12);
        xlabel('TMS timing','fontsize',12);
        title('Control: mised effect model estimate','fontsize',14);
        
        % simple means
        figure;
        lineplot(DP.(timing_var),DP.perc,'split',DP.week,'subset',~isnan(DP.perc),...
            'leg',Weeks,'markersize',8,'markercolor',color_w,'markerfill',color_w,...
            'linecolor',color_w,'errorcolor',color_w,'linewidth',1.5,'errorcap',0,'errorwidth',1.5);
        drawline(0,'dir','horz');
        set(gca,'ylim',[-0.6 0.4]);
        ylabel('(CT-UT)/UT','fontsize',12);
        title('Patient: simple mean','fontsize',14);
        if strcmp(timing_type,'real')
            set(gca,'XTickLabel',xtlabel);
        end
        
        figure;
        %lineplot(DC.tms_schedule,DC.perc,'split',DC.week,'subset',~isnan(DC.perc),...
        lineplot(DC.(timing_var),DC.perc,'subset',~isnan(DC.perc),...
            'leg',Weeks,'markersize',8,'markercolor',color_w,'markerfill',color_w,...
            'linecolor',color_w,'errorcolor',color_w,'linewidth',1.5,'errorcap',0,'errorwidth',1.5);
        drawline(0,'dir','horz');
        set(gca,'ylim',[-0.6 0.4]);
        ylabel('(CT-UT)/UT','fontsize',12);
        xlabel('TMS timing','fontsize',12);
        title('Control: simple mean','fontsize',14);
        if strcmp(timing_type,'real')
            set(gca,'XTickLabel',xtlabel);
        end
        
    case 'plot_dIHI'
        interval=varargin{1};  % RT interval to look at, e.g. [0.2 0.95] is the interval between 0.20 & 0.95%RT
                               % [0.35,0.88] is the mean(0.2,0.5) vs. mean(0.8,0.95)
        D=load(fullfile(analysisDir,'summary_tms_time.mat'));
        [XP,R]=pivottable(D.tms_schedule,D.week,D.perc,'nanmean','subset',D.control==0);
        XP(5,:) = mean(XP(1:2,:)); % mean over 2 early tms timings
        XP(6,:) = mean(XP(3:4,:)); % mean over 2 late tms timings
        R(5) = 0.35;
        R(6) = 0.88;
        
        [XC,R]=pivottable(D.tms_schedule,D.week,D.perc,'nanmean','subset',D.control==1);
        XC(5,:) = mean(XC(1:2,:)); % mean over 2 early tms timings
        XC(6,:) = mean(XC(3:4,:)); % mean over 2 late tms timings
        
        figure;
        color=[[0 0 0];[0.5 0 0];[1 0 0];[1 0.5 0];[1 1 0]];
        indx=(1:5);
        bh=bar([X(find(R==interval(2)),1)-X(find(R==interval(1)),1);...
                X(find(R==interval(2)),2)-X(find(R==interval(1)),2);...
                X(find(R==interval(2)),3)-X(find(R==interval(1)),3);...
                X(find(R==interval(2)),4)-X(find(R==interval(1)),4);...
                X(find(R==interval(2)),5)-X(find(R==interval(1)),5)],'EdgeColor','none');
        set(gca,'XTickLabel',weeks);
        bhch=get(bh,'Children');
        set(bhch,'CData',indx);
        colormap(color);
        ylabel('Delta IHI (88-35%RT)','FontSize',12);
        set(gca,'ylim',[0 0.5]);
        if control
            title('Control','Fontsize',14);
        else
            title('Patient','Fontsize',14);
        end
        
    case 'mixed_model'               % Mixed model for mean across different TMS timing across weeks
        D     = varargin{1}; % all data
        varin = varargin{2}; % variable to be estimated, vector of values
        
        fig = 0;
        ff_week = 0;
        vararginoptions(varargin(3:end),{'fig','ff_week'});
        T = getrow(D,~isnan(varin));     % Find the observed data
        v = varin(~isnan(varin));
        Subj = pivottable(T.SubjN,T.control,T.SubjN,'mean');
        
        % Build fixed effects and random effects design matrix
        if ff_week
            M   = mixed_model('make_repeatedMeasures',T.SubjN,{T.week,'reduced',0});
            M.S = Subj;
            fflevels = 0;
            R.DM = weeks;
        else
            fflevels = unique(T.splitvar);        
            if length(fflevels)==2
                M = mixed_model('make_repeatedMeasures',T.SubjN,{T.splitvar,'reduced',1},...
                    {T.week.*(T.splitvar==fflevels(1)),'reduced_p',0},... % Recovery of non-paretic (or left)
                    {T.week.*(T.splitvar==fflevels(2)),'reduced_p',0});   % Recovery of paretic (or right)
            elseif length(fflevels)==3
                M = mixed_model('make_repeatedMeasures',T.SubjN,{T.splitvar,'reduced',1},...
                    {T.week.*(T.splitvar==fflevels(1)),'reduced_p',0},...  % Recovery of non-paretic (or left)
                    {T.week.*(T.splitvar==fflevels(2)),'reduced_p',0},...  % Recovery of paretic (or right)
                    {T.week.*(T.splitvar==fflevels(3)),'reduced_p',0});    % Recovery of paretic hand in severe patients
            elseif length(fflevels)==4
                M = mixed_model('make_repeatedMeasures',T.SubjN,{T.splitvar,'reduced',1},...
                    {T.week.*(T.splitvar==fflevels(1)),'reduced_p',0},...  
                    {T.week.*(T.splitvar==fflevels(2)),'reduced_p',0},...  
                    {T.week.*(T.splitvar==fflevels(3)),'reduced_p',0},...  
                    {T.week.*(T.splitvar==fflevels(4)),'reduced_p',0});    
            elseif length(fflevels)==5
                M = mixed_model('make_repeatedMeasures',T.SubjN,{T.splitvar,'reduced',1},...
                    {T.week.*(T.splitvar==fflevels(1)),'reduced_p',0},...  
                    {T.week.*(T.splitvar==fflevels(2)),'reduced_p',0},...  
                    {T.week.*(T.splitvar==fflevels(3)),'reduced_p',0},...  
                    {T.week.*(T.splitvar==fflevels(4)),'reduced_p',0},...  
                    {T.week.*(T.splitvar==fflevels(5)),'reduced_p',0});    
            end
        end
        M   = mixed_model('estimate_reml',M,v);
        M.S = Subj;
        
        % Now get the mean response that we would like to investigate -
        % give design matrix
        if length(fflevels) > 1
            R.splitF = [];
            for i=1:length(fflevels)
                R.splitF = [R.splitF; ones(5,1)*fflevels(i)];
                R.week   = repmat(weeks,length(fflevels),1);
            end
            R.DM = R.splitF;
            for i=1:length(fflevels)
                R.DM = [R.DM R.week.*(R.splitF==fflevels(i))];
            end
        end
        [R.y,R.yvar,V,X] = mixed_model('meanResponse',R.DM,M);
                
        if (fig)
            figure(1);
            subplot(1,2,1);
            lineplot(T.week,varin,'split',T.splitvar,...
                'markersize',10,'markerfill',paretic_color,'linewidth',1.5);
            title('Simple mean estimate','FontSize',16);
            xlabel('Weeks','FontSize',12);
            
            subplot(1,2,2);
            lineplot(R.week,R.y,'split',R.splitF,'errorval',sqrt(R.yvar),...
                'markersize',10,'markerfill',paretic_color,'linewidth',1.5);
            title('Random effects mean','FontSize',16);
            xlabel('Weeks','FontSize',12);
            set(gcf,'PaperPositionMode','auto');
                        
            keyboard;
        end;
        varargout={R,M};

    % OLD cases
end