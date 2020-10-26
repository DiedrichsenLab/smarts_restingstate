function clean_ihi_data(ID,Week)
% check "IHI_list.xlsx" for details

data_dir = '/Volumes/dayu/McDonnell/TMS/data/preprocess/IHI/';
subj_dir = fullfile(data_dir,ID,Week);

% Block 6 was run using B1 protocol, copy experimental data to B6
if strcmp(ID,'JHU2650')&&strcmp(Week,'W52')
    p = load('IHI_JHU2650_W52_paradigm');
    p.tms.states(121:end)    = p.tms.states(1:24);
    p.tms.tms_time(121:end)  = p.tms.tms_time(1:24);
    p.tms.scheduled(121:end) = p.tms.scheduled(1:24);
    save('nIHI_JHU2650_W52_paradigm.mat','-struct','p');
end

% break IHIP1 to two blocks
if strcmp(ID,'UZ3238')&&strcmp(Week,'W4')
    cd(subj_dir);
    meps = load('IHIP1.UZ3238.FDIR.W4.MEP.mat');
    names = fieldnames(meps.Summary.channel);
    for f=1:length(names)
        for ch=1:2
            mep1.Summary.channel(ch).(names{f}) = meps.Summary.channel(ch).(names{f})(1:24);
            mep2.Summary.channel(ch).(names{f}) = meps.Summary.channel(ch).(names{f})(25:48);
        end
    end
    names = fieldnames(meps.MEPs.channel);
    for f=1:length(names)
        for ch=1:2
            mep1.MEPs.channel(ch).(names{f}) = meps.MEPs.channel(ch).(names{f})(1:24);
            mep2.MEPs.channel(ch).(names{f}) = meps.MEPs.channel(ch).(names{f})(25:48);
        end
    end
    waves = load('IHIP1.UZ3238.FDIR.W4.mat');
    name1 = fieldnames(waves);
    name2 = fieldnames(waves.(name1{1}));
    for f=1:length(name2)
        wave1.(name1{1}).(name2{f}) = waves.(name1{1}).(name2{f});
        wave2.(name1{1}).(name2{f}) = waves.(name1{1}).(name2{f});
    end
    wave1.(name1{1}).frames    = 24;
    wave2.(name1{1}).frames    = 24;
    wave1.(name1{1}).frameinfo = waves.(name1{1}).frameinfo(1:24);
    wave2.(name1{1}).frameinfo = waves.(name1{1}).frameinfo(25:48);
    wave1.(name1{1}).values    = waves.(name1{1}).values(:,:,1:24);
    wave2.(name1{1}).values    = waves.(name1{1}).values(:,:,25:48);
    save('nIHIP1.UZ3238.FDIR.W4.MEP.mat','-struct','mep1');
    save('nIHIP2.UZ3238.FDIR.W4.MEP.mat','-struct','mep2');
    save('nIHIP1.UZ3238.FDIR.W4.mat','-struct','wave1');
    save('nIHIP2.UZ3238.FDIR.W4.mat','-struct','wave2');
end

% Block 6 contains only 17 trials, remove the last 7 trials from paradigm file
if strcmp(ID,'UZ3239')&&strcmp(Week,'W12')
    cd(subj_dir);
    p = load('IHI_UZ3239_W12_paradigm');
    p.tms.TN(end-6:end)         = [];
    p.tms.go_onset(end-6:end)   = [];
    p.tms.states(end-6:end)     = [];
    p.tms.tms_time(end-6:end)   = [];
    p.tms.scheduled(end-6:end)  = [];
    p.tms.ITI(end-6:end)        = [];
    p.exp_params.ITIs(end-6:end)= [];
    save('nIHI_UZ3239_W12_paradigm.mat','-struct','p');
end    

% Block 5 used B1 protocol, correct tms_time in paradigm file, copy b1->b5
if strcmp(ID,'UZ3246')&&strcmp(Week,'W12')
    cd(subj_dir);
    p = load('IHI_UZ3246_W12_paradigm');
    p.tms.tms_time(97:120)   = p.tms.tms_time(1:24);
    p.tms.scheduled(97:120)  = p.tms.scheduled(1:24);
    p.tms.ITI(97:120)        = p.tms.ITI(1:24);
    p.exp_params.ITIs(97:120)= p.exp_params.ITIs(1:24);
    save('nIHI_UZ3246_W12_paradigm.mat','-struct','p');    
end

% Block 1 contains one extra trial at the end, remove from .mat & MEP files
if strcmp(ID,'UZP1007')&&strcmp(Week,'W4')
    cd(subj_dir);
    meps = load('IHIP1.UZP1007.FDIL.W4.MEP.mat');
    names = fieldnames(meps.Summary.channel);
    for f=1:length(names)
        for ch=1:2
            meps.Summary.channel(ch).(names{f})(end) = [];
        end
    end
    waves = load('IHIP1.UZP1007.FDIL.W4.mat');
    name = fieldnames(waves);
    waves.(name{1}).frames    = 24;
    waves.(name{1}).frameinfo(end) = [];
    waves.(name{1}).values(:,:,end)= [];
    save('nIHIP1.UZP1007.FDIL.W4.MEP.mat','-struct','meps');
    save('nIHIP1.UZP1007.FDIL.W4.mat','-struct','waves');
end

% correct for old script in UZ data
if (strcmp(ID,'UZ3166')&&strcmp(Week,'W24'))||...
   (strcmp(ID,'UZP1001')&&strcmp(Week,'W1'))||...
   (strcmp(ID,'UZP1001')&&strcmp(Week,'W4'))||...
   (strcmp(ID,'UZP1002')&&strcmp(Week,'W52'))||...
   (strcmp(ID,'UZP1004')&&strcmp(Week,'W1'))
    cd(subj_dir);
    correct_oldscript(ID,Week);
end

if (strcmp(ID,'UZP1002')&&strcmp(Week,'W4'))
    % frist remove block 3 data from paradigm file
    p = load('IHI_UZP1002_W4_paradigm.mat');
    
    % then correct for old script
    correct_oldscript(ID,Week);
end

function correct_oldscript1(ID,Week)
    p = load(['IHI_' ID '_' Week '_paradigm']);
    blocks = dir(['IHI*' Week '.mat']); % assume blocks all named & ordered correctly
    muscle = blocks(1).name(length(ID)+8:length(ID)+11);
    ntrials = 24;
    for b=1:length(blocks)
        filename = ['IHIP' num2str(b) '.' ID '.' muscle '.' Week '.mat'];
        waves  = load(filename);
        name   = fieldnames(waves);
        states = [];
        for t=1:ntrials
            states(t,1) = waves.(name{1}).frameinfo(t).state;
        end
        % remove mis-matched trials
        trials = ntrials*(b-1);
        idx = find(p.tms.states((trials+1):(trials+24))~=states);
        p.tms.states(trials+idx)   = NaN;
        p.tms.tms_time(trials+idx) = NaN;
        p.tms.scheduled(trials+idx)= NaN;
        waves.(name{1}).frames = waves.(name{1}).frames-length(idx);
        waves.(name{1}).frameinfo(idx) = [];
        waves.(name{1}).values(:,:,idx)= [];
        save(['n' filename],'-struct','waves');
    end
    keyboard;
    idx=find(isnan(p.tms.states));
    p.tms.states(idx)   = [];
    p.tms.tms_time(idx) = [];
    p.tms.scheduled(idx)= [];
    save(['nIHI_' ID '_' Week '_paradigm.mat'],'-struct','p');

function correct_oldscript(ID,Week)
    p = load(['IHI_' ID '_' Week '_paradigm']);
    blocks = dir(['IHIP*' Week '.mat']); % assume blocks all named & ordered correctly
    muscle = blocks(1).name(length(ID)+8:length(ID)+11);
    ntrials = 24;
    for b=1:length(blocks)
        meps = load(['IHIP' num2str(b) '.' ID '.' muscle '.' Week '.MEP.mat']);
        trials = ntrials*(b-1);
        % remove mis-matched trials
        idx = find(p.tms.states((trials+1):(trials+24))~=meps.Summary.channel(1).state);
        p.tms.states(trials+idx)   = meps.Summary.channel(1).state(idx);
        % corret RT (+8ms), because used EMG delay 320ms, which should be 312ms
        
    end
    keyboard;
    save(['nIHI_' ID '_' Week '_paradigm.mat'],'-struct','p');
