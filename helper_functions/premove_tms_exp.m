function premove_tms_exp
% Premovement IHI/SICI experiment for McDonnell project
% Instructions for setting up in "premovement_TMS_setup.docx"
% (c) Jing Xu, jing.xu@jhmi.edu
% 07/27/2012
% 03/04/2014 Jing Xu - added 0-2sec jittering to ITI

Screen('CloseAll');
clear all;
addpath(genpath(pwd));
global gExp;
base_dir = pwd;
data_dir = fullfile(base_dir,'data');

% get subject information: subject number
accepted = false;
while (~accepted)
    fprintf('Welcome to TMS experiment.\n');
    D.subj = input('Enter subject ID: ','s');
    if isempty(ls(sprintf('%s\\*%s*.mat',data_dir,D.subj)))
        accepted = true;
    else
        answer = input('\nSubject already exist. Continue? (Y/N)','s');
        if strcmpi(answer,'Y')
            accepted = true;
        end
    end
end
rt_data_file = fullfile(data_dir, ['RT_' D.subj '_paradigm.mat']);
n_rttest = 30;
n_protocols = 18;  % number of trials per type of tms pulse
n_blocks = 6; % 3 - 48 trials per block; 6 - 24 trials per block
n_states = 2;
%gExp.ITI = 5;
stim_duration = 2; % time to present go cue, in sec
tms_timing = [0.2 0.5 0.8 0.95];
ITI_jitter = linspace(0,2,n_blocks*24);
ITIs = 5+ITI_jitter(randperm(n_blocks*24));
ITI_jitter = linspace(0,2,n_rttest);
ITIs_rt = 5+ITI_jitter(randperm(n_rttest));
gExp.pause = false;
% for saving data: subject info, go stim onset time, tms timing and trial type
D.rt.TN = (1:n_rttest)';
D.rt.go_onset = zeros(n_rttest,1);
n_timing = length(tms_timing);
n_tms = n_protocols*n_states*n_timing;
D.tms.TN = (1:n_tms)';
D.tms.go_onset = zeros(n_tms,1);
D.tms.states = zeros(n_tms,1);
D.tms.tms_time = zeros(n_tms,1);
D.tms.scheduled = zeros(n_tms,1);
D.exp_params.n_rttest = n_rttest;  % book-keeping experimental params
D.exp_params.n_protocols = n_protocols;
D.exp_params.n_blocks = n_blocks;
D.exp_params.n_states = n_states;
%D.exp_params.ITI = gExp.ITI;
D.exp_params.ITIs = ITIs;
D.exp_params.stim_duration = stim_duration;
D.exp_params.tms_timing = tms_timing;

init_graphics; % prepare graphics & present the welcome message
isdone = false;
while(~isdone)
    gExp.pause = false;
    str=input('\nEXP>','s');
    [command,param] = strtok(str);
    open_parallel;  % opens parallel port
    switch command
        case 'test_rt'  % Part I: 30 trials of simple Go RT, trigger EMG at Go onset
            update_graphics('test_rt', false);
            fprintf('\nPress space bar to start; press q to break.');
            FlushEvents('keyDown');
            KbStrokeWait;
            for i = 1:n_rttest
                avail = CharAvail;
                if avail > 0
                    char = GetChar;
                    if any(char == 'q')
                        fprintf('\nUser broke out of experiment.\n');
                        gExp.pause=true;
                        break;
                    end
                end
                ITI = ITIs_rt(i);
                fprintf(['ITI: ' num2str(ITI)]);
                [D.rt.go_onset(i)] = run_trial(stim_duration,0,ITI);
            end
            save(rt_data_file,'-struct','D');
            update_graphics('rt_done',false);
        case 'set'                          % Sets any global variable the value specified
            [varName, value]= strtok(param);
            if isfield(gExp,varName)
                gExp.(varName)= str2num(value);
            else
                fprintf('%s is not a valid variable name.',varName);
            end
        case 'test_ts'  % 50% RT for testing TS intensities
            test_ts_isdone = false;
            update_graphics('tms',false);
            go_rt = input('\nEnter EMG RT: ');
            c = 1;
            while ~test_ts_isdone
                update_graphics('break', false);
                FlushEvents('keyDown');
                fprintf('\nPress any key to continue; press q to quit.\n');
                KbStrokeWait;
                for i = 1:10 % 10 trials per block
                    avail = CharAvail;
                    if avail > 0
                        char = GetChar;
                        if any(char == 'q')
                            fprintf('\nUser broke out of experiment.\n');
                            gExp.pause=true;
                            test_ts_isdone = true;
                            break;
                        end                            
                    end
                    if c>length(ITIs_rt)
                        c = c - length(ITIs_rt);
                    end
                    ITI = ITIs_rt(c); c=c+1;
                    fprintf(['ITI: ' num2str(ITI)]);
                    run_trial(stim_duration,0.5*go_rt,ITI);
                end
            end
        case 'tms'  % pre-movement tms, for IHI/SICI/LICI
        % Part II: 144 trials of TMS: 72 single-pulse & 72 paired-pulse
        % show screen for reminders: Signal, EMG, TMS intensity
        % show screen for instructions
            D.session = input('Enter session (IHI/SICI): ','s');
            data_file = fullfile(data_dir, [D.session '_' D.subj '_paradigm.mat']);
            update_graphics('tms',false);
            % load sequence protocols for Signal
            load(fullfile('protocols', 'protocols.mat'));
            go_rt = input('\nEnter EMG RT: ');
            D.tms.pre_gort = go_rt;
            % start trial loop: 
            b = 1;
            while (b<=n_blocks)
                gExp.pause=false;
                update_graphics('break', false);
                fprintf('\nBlock %d\n',b);
                FlushEvents('keyDown');
                fprintf('\nPress any key to continue; press q to quit.\n');
                KbStrokeWait;
                avail = CharAvail;
                if avail > 0
                    char = GetChar;
                    if any(char == 'q')
                        fprintf('\nUser broke out of experiment.\n');
                        gExp.pause=true;
                        b = b-1;
                        break;
                    end
                end
                for i = 1:(n_protocols/n_blocks)
                    avail = CharAvail;
                    if avail > 0
                        char = GetChar;
                        if any(char == 'q')
                            fprintf('\nUser broke out of experiment.\n');
                            gExp.pause=true;
                            b = b-1;
                            break;
                        end
                    end
                    single_timing = tms_timing(randperm(n_timing));
                    pair_timing = tms_timing(randperm(n_timing));
                    st = 1; pt = 1; % pointers for single-pulse and paired-pulse timing
                    for j = 1:(n_timing*n_states)  % number of trials per protocol
                        avail = CharAvail;
                        if avail > 0
                            char = GetChar;
                            if any(char == 'q')
                                fprintf('\nUser broke out of experiment.\n');
                                gExp.pause=true;
                                b = b-1;
                                break;
                            end
                        end
                        p = (b-1)*(n_protocols/n_blocks)+i; % protocol number
                        t = n_timing*n_states*(p-1)+j; % trial number
                        state = protocols(p,j);
                        ITI = ITIs(t);
                        D.tms.states(t) = state;
                        D.tms.ITI(t) = ITI;
                        if state == 1
                            fprintf(['\nstate: ' num2str(state)]);
                            fprintf(['  ITI: ' num2str(ITI)]);
                            fprintf(['  timing: ' num2str(single_timing(st)) ' ']);
                            D.tms.scheduled(t) = single_timing(st);
                            [D.tms.go_onset(t),D.tms.tms_time(t)] = ...
                                run_trial(stim_duration,single_timing(st)*go_rt,ITI);
                            st = st + 1;
                        else
                            fprintf(['\nstate: ' num2str(state)]);
                            fprintf(['  ITI: ' num2str(ITI)]);
                            fprintf(['  timing: ' num2str(pair_timing(pt)) ' ']);
                            D.tms.scheduled(t) = pair_timing(pt);
                            [D.tms.go_onset(t),D.tms.tms_time(t)] = ...
                                run_trial(stim_duration,pair_timing(pt)*go_rt,ITI);
                            pt = pt + 1;
                        end
                    end
                    if gExp.pause  % if user broke out of the exp, back to >EXP environment
                       break;
                    end                     
                end
               b = b+1;
            end
            save(data_file,'-struct','D');
        case {'quit','exit'}
            isdone=true;
            update_graphics('end_exp',false);
            WaitSecs(2);
        case ''                             % just hit return: do nothing
        otherwise                           % Unknown command: give error message
            fprintf('Unknown command\n');
    end
end
Screen('CloseAll');
clear mex;
clear('gExp');

function [stim_onset,trigger_time] = run_trial(stim_duration, tms_time, ITI)
%global gExp;
% If tms_time = 0, trigger at go stimulus onset time
update_graphics('running',false);
WaitSecs(ITI);

% present Go stim
stim_onset = update_graphics('running',true);
WaitSecs(tms_time/1000);
trigger_time = trigger_signal;  % get actual trigger time
WaitSecs(stim_duration-(tms_time/1000));
% stop Go stim after 1 sec
stim_offset = update_graphics('running',false);
fprintf('\n\nstim time: %f,\ntms scheduled: %f,\nactual tms: %f\n',...
    (stim_offset-stim_onset), tms_time, trigger_time-stim_onset);

function varargout = open_parallel
global gExp;
[gExp.portData, gExp.portStatus, gExp.portControl] = ptb_init_parallel_port('Open');
putvalue(gExp.portData,0);

if nargout == 1
    varargout = {gExp.portData};
elseif nargout == 2
    varargout = {gExp.portData,gExp.portStatus};
elseif nargout == 3
    varargout = {gExp.portData,gExp.portStatus,gExp.portControl};
end

function t = trigger_signal
global gExp;
putvalue(gExp.portData(1),1);
t = GetSecs;
WaitSecs(0.001);
putvalue(gExp.portData(1),0);

function init_graphics
global gExp;
gExp.myColor=[0,0,0;...   % 1: Black
    255,255,255;...		  % 2: White
    0,200,0;...   		  % 3: Green
    150,0,0;...			  % 4: Red
    50,50,50;...	      % 5: gray
    30,30,30;...          % 6: darkgray
    0,0,200];             % 7: blue

% openning up a blank black screen
Screen('Preference', 'VisualDebugLevel', 0);
n_monitors = size(Screen('Screens'),2);
if (n_monitors == 3)
    [gExp.wPtr,gExp.ScreenRect]=Screen('OpenWindow',2,gExp.myColor(1,:)); % projecting on second screen
    scpixels = [1 1 gExp.ScreenRect(3:4)];
    monitor_size = scpixels / get(0,'ScreenPixelsPerInch') * 0.0254;  % convert to meters
    [gExp.ScreenWidth, gExp.ScreenHeight]=Screen('WindowSize', gExp.wPtr);
else
    % get monitor width and height
    scpixels = get(0,'ScreenSize');
    monitor_size = scpixels / get(0,'ScreenPixelsPerInch') * 0.0254;  % convert to meters
    gExp.ScreenWidth = 1000;  % in pixels
    gExp.ScreenHeight =1000;
    x_offset = 10;
    y_offset = 10;
    [gExp.wPtr,gExp.ScreenRect]=Screen('OpenWindow',0,gExp.myColor(1,:),...
        [x_offset y_offset gExp.ScreenWidth gExp.ScreenHeight]);
end
xorg = gExp.ScreenWidth/ 2;
yorg = gExp.ScreenHeight / 2;
pix_meter_x = scpixels(3) / monitor_size(3);  % used for conversion from meter to pixels
pix_meter_y = scpixels(4) / monitor_size(4);

% these can be modified to other shapes later
gExp.gocue.type = 'circle';
gExp.gocue.x = 0;  % meters away from center
gExp.gocue.y = 0;  % meters away from center
gExp.gocue.radius = 0.01;  % in meters
gExp.gocue.color = gExp.myColor(3,:);
tlx = gExp.gocue.x - gExp.gocue.radius;
tly = gExp.gocue.y - gExp.gocue.radius;
brx = gExp.gocue.x + gExp.gocue.radius;
bry = gExp.gocue.y + gExp.gocue.radius;
gExp.gocue.outrect(1) = floor(xorg + (tlx*pix_meter_x));
gExp.gocue.outrect(2) = floor(yorg + (tly*pix_meter_y));
gExp.gocue.outrect(3) = floor(xorg + (brx*pix_meter_x));
gExp.gocue.outrect(4) = floor(yorg + (bry*pix_meter_y));

gExp.text.font = 'Arial';
gExp.text.size = 48;
gExp.text.style = 0;
gExp.text.message{1} = 'Welcome to the TMS studies!';
gExp.text.message{2} = 'Starting Pre-Test.';
gExp.text.message{3} = 'Pre-Test done.';
gExp.text.message{4} = 'Starting TMS test.';
gExp.text.message{5} = 'Get ready for the next block.';
gExp.text.message{6} = 'Thank you for participating!';
Screen('TextFont',gExp.wPtr,gExp.text.font);
Screen('TextSize',gExp.wPtr,gExp.text.size);
Screen('TextStyle',gExp.wPtr,gExp.text.style);
bounds=Screen('TextBounds',gExp.wPtr,gExp.text.message{1});
out_rect = CenterRectOnPoint(bounds,0.5*gExp.ScreenWidth,0.5*gExp.ScreenHeight);
Screen('DrawText',gExp.wPtr,gExp.text.message{1},out_rect(1),out_rect(2),gExp.myColor(2,:));
Screen('Flip', gExp.wPtr);


function t = update_graphics(exp_state,stim_on)
global gExp;
switch exp_state
    case 'test_rt'
        Screen('TextFont',gExp.wPtr,gExp.text.font);
        Screen('TextSize',gExp.wPtr,gExp.text.size);
        Screen('TextStyle',gExp.wPtr,gExp.text.style);
        bounds=Screen('TextBounds',gExp.wPtr,gExp.text.message{2});
        out_rect = CenterRectOnPoint(bounds,0.5*gExp.ScreenWidth,0.5*gExp.ScreenHeight);
        Screen('DrawText',gExp.wPtr,gExp.text.message{2},out_rect(1),out_rect(2),gExp.myColor(2,:));
        [t1,t2,t] = Screen('Flip', gExp.wPtr);
    case 'rt_done'
        Screen('TextFont',gExp.wPtr,gExp.text.font);
        Screen('TextSize',gExp.wPtr,gExp.text.size);
        Screen('TextStyle',gExp.wPtr,gExp.text.style);
        bounds=Screen('TextBounds',gExp.wPtr,gExp.text.message{3});
        out_rect = CenterRectOnPoint(bounds,0.5*gExp.ScreenWidth,0.5*gExp.ScreenHeight);
        Screen('DrawText',gExp.wPtr,gExp.text.message{3},out_rect(1),out_rect(2),gExp.myColor(2,:));
        [t1,t2,t] = Screen('Flip', gExp.wPtr);
    case 'tms'
        Screen('TextFont',gExp.wPtr,gExp.text.font);
        Screen('TextSize',gExp.wPtr,gExp.text.size);
        Screen('TextStyle',gExp.wPtr,gExp.text.style);
        bounds=Screen('TextBounds',gExp.wPtr,gExp.text.message{4});
        out_rect = CenterRectOnPoint(bounds,0.5*gExp.ScreenWidth,0.5*gExp.ScreenHeight);
        Screen('DrawText',gExp.wPtr,gExp.text.message{4},out_rect(1),out_rect(2),gExp.myColor(2,:));
        [t1,t2,t] = Screen('Flip', gExp.wPtr);
    case 'running'
        if stim_on
            Screen('FillOval',gExp.wPtr,gExp.gocue.color,gExp.gocue.outrect);
            [t1,t2,t] = Screen('Flip',gExp.wPtr,0,0,2);
        else
            % clear up the screen
            Screen('FillRect', gExp.wPtr, 0, gExp.ScreenRect);
            %Screen('FillOval',gExp.wPtr,[0 0 0],gExp.gocue.outrect);
            [t1,t2,t] = Screen('Flip', gExp.wPtr);
        end
    case 'break'
        Screen('TextFont',gExp.wPtr,gExp.text.font);
        Screen('TextSize',gExp.wPtr,gExp.text.size);
        Screen('TextStyle',gExp.wPtr,gExp.text.style);
        bounds=Screen('TextBounds',gExp.wPtr,gExp.text.message{5});
        out_rect = CenterRectOnPoint(bounds,0.5*gExp.ScreenWidth,0.5*gExp.ScreenHeight);
        Screen('DrawText',gExp.wPtr,gExp.text.message{5},out_rect(1),out_rect(2),gExp.myColor(2,:));
        [t1,t2,t] = Screen('Flip', gExp.wPtr);
    otherwise  % end of experiment
        Screen('TextFont',gExp.wPtr,gExp.text.font);
        Screen('TextSize',gExp.wPtr,gExp.text.size);
        Screen('TextStyle',gExp.wPtr,gExp.text.style);
        bounds=Screen('TextBounds',gExp.wPtr,gExp.text.message{6});
        out_rect = CenterRectOnPoint(bounds,0.5*gExp.ScreenWidth,0.5*gExp.ScreenHeight);
        Screen('DrawText',gExp.wPtr,gExp.text.message{6},out_rect(1),out_rect(2),gExp.myColor(2,:));
        [t1,t2,t] = Screen('Flip', gExp.wPtr);
end
