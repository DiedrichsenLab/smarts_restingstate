function [Summary,MEPs,optflags] = process_mep2(dd,flags);

% USAGE: [Summary,MEPs,flags] = process_mep(DataIn,flags);
%
% DataIn can be either a filename (full path), a datastructure
% (loaded from a file exported from Signal 4.02), empty, or not
% supplied.  If empty or not supplied, user will be prompted to
% select a *.mat file (exported from Signal 4.02).
%
% flags is a structure for analysis flags.  If not supplied,
% defaults will be used.  If any fields are missing, defaults will
% be used.  Fields as follows:
%
%   channel -- data channel to analyze.  Default is -1, which
%      prompts the user to specify manually.
%   gain -- gain of Viking system.  Default is -1, which prompts
%      the user to specify manually.
%   background_noise -- boolean, to measure background noise.
%      Default is 0.
%   handfilter -- boolean, to allow user to filter through MEP
%      traces.  Default is 1 (highly recommended!!!)
%
% Summary will be the standard dsave/dload file.
% MEPs will have all the same data, as well as MEP traces.
% flags as output will store a copy of the option flags used.
 
% Essentially a streamlining of mep_any, by Joe Galea
%   -- John Schlerf, August 2010


% set up default option flags:
    
% optflags = struct('channel',-1,'gain',-1,'background_noise',1,...
%                   'mepTime',[1:700],'preTime',[51:200],...
%                   'totalTime',[1:750],'handfilter',1);
%               
% optflags = struct('channel',-1,'gain',-1,'background_noise',1,...
%                   'mepTime',1250+[50:995],'preTime',[51:200],...
%                   'totalTime',[1:3750],'handfilter',1);
              
 optflags = struct('channel',-1,'gain',-1,'background_noise',1,...
                   'mepTime',[2600:3500],'preTime',[2600:2750],...
                   'totalTime',[1:5000],'handfilter',1);


% parse my input:
if exist('dd','var')
    if ~exist('flags','var') 
        if isstruct(dd) && any(ismember(fields(dd),fields(optflags)))
            flags = dd;
            dd = [];
        end
    end
end

if exist('flags','var')
    f = fields(flags);
    for i=1:length(f)
        optflags = setfield(optflags,f{i},getfield(flags,f{i}));
    end
end

if ~exist('dd','var') || isempty(dd)
    disp('Load data');
    [file] = find_file('*.mat*'); 
elseif exist(dd,'file')
    file = dd;
elseif ischar(dd)
    disp('Cannot find specified file, please reselect.');
    [file] = find_file('*.mat*');
end


% load if needed:
if exist('file','var')
    
    dd = load(file); 
    f = fields(dd);
    dd = getfield(dd,f{1});
    
end

n = dd.frames;
channels=dd.chans;
good_trial=ones(1,n);
thresh = 0.25;
temp = zeros(2,4);
data_channel = optflags.channel;
if data_channel == -1
    data_channel = input('please input the data channel: ');
end
gain = optflags.gain;
if gain == -1
    gain = input('please input the gain: ');
end

background_noise = optflags.background_noise;
mepTime = optflags.mepTime;
totalTime = optflags.totalTime;
preTime = optflags.preTime;

%trim trial to onset of pulse + 200 msecs.
for i=1:n
    [pulse(i),time1(i)]=max(dd.values(:,1,i));  
end

avg_pulse=mean(pulse);
avg_pulse_time=mean(time1);

%% process all data that I want, trial by trial

fighandle = figure; 
mepAxes = axes('position',[0.05,0.55,0.9,0.4]);
hold on;
if background_noise
    wholeTraceAxes = axes('position',[0.05,0.05,0.6,0.4]);
    hold on;
    bgAxes = axes('position',[0.75,0.05,0.2,0.4]);
    hold on;
else
    wholeTraceAxes = axes('position',[0.05,0.05,0.9,0.4]);
    hold on;
    bgAxes = 0;
end


Summary = struct;
MEPs = struct;

frameinfoFields = fields(dd.frameinfo(1));
i=1;

while i <= n
    
    skipSwitch = 0;
    
    M = struct;
    M.trace = {(gain/200)*squeeze(dd.values(mepTime,data_channel,i))};
    M.traceWhole = {(gain/200)*squeeze(dd.values(totalTime,data_channel,i))};
    M.time = {mepTime};
    % Guess min and max of each trial
    [D.peak,D.peakTime] = max(M.trace{1});
    [D.trough,D.troughTime] = min(M.trace{1});
    
    % Pass frameinfo values into output:
    for f=1:length(frameinfoFields)
        thisValue = getfield(dd.frameinfo(i),frameinfoFields{f});
        D = setfield(D,frameinfoFields{f},thisValue);
        M = setfield(M,frameinfoFields{f},thisValue);
    end
    
    if background_noise
        D.rms = sqrt(mean(dd.values(preTime,data_channel,i).^2));
        D.pre_act_trial_avg = mean(abs(dd.values(preTime, ...
                                                data_channel,i)));
        D.pre_act_trial_std = std(abs(dd.values(preTime, ...
                                                data_channel,i)));
        M.bg = {squeeze(dd.values(preTime,data_channel,i))};
        M.bg_rec = {abs(M.bg{1})};
    end
    
    axes(wholeTraceAxes);
    cla
    plot(M.traceWhole{1});
    plot(D.peakTime+mepTime(1),D.peak,'R*');
    plot(D.troughTime+mepTime(1),D.trough,'G*');
    
%     if bgAxes
%         axes(bgAxes)
%         cla
%         bar(1,D.pre_act_trial_avg,'b');
%         errorbars(1,D.pre_act_trial_avg,D.pre_act_trial_std);
%         bar(2,D.peak-D.trough,'r');
%     end    
    
    axes(mepAxes);
    % plot MEP, and adjust max/min if necessary:
    cla
    plot(M.trace{1});
    plot(D.peakTime,D.peak,'R*');
    plot(D.troughTime,D.trough,'G*');
    
    %[x,y,button]=ginput(1);

    if optflags.handfilter
        while 1, %allows for 3 button mice now
        
            try
                [x,y,button]=ginput(1);
                x=round(x);
            catch
                button=2;
                lasterr('');        
            end
            
            switch button
                
              case {1}
                %keyboard
                %D.peakTime=x;
                % Use the nearest maximal value, not the mouse value:
                dt = [0; sign(diff(M.trace{1}))];
                zCrossLoc=find(diff(dt));
                distSq = (zCrossLoc - x).^2;
                [minDist,which] = min(distSq);
                D.peakTime = zCrossLoc(which);
                D.peak=M.trace{1}(zCrossLoc(which));
              case (2)
                D.isGood = ...
                    isequal(questdlg('Is this a good trial?',...
                                     'good trial','yes','no','yes'),'yes');
                break;    
              case {3}
                %keyboard
                %D.troughTime=x;
                % Use the nearest minimum value, not the mouse value:
                dt = [0; sign(diff(M.trace{1}))];
                zCrossLoc=find(diff(dt));
                distSq = (zCrossLoc - x).^2;
                [minDist,which] = min(distSq);
                D.troughTime = zCrossLoc(which);
                D.trough=M.trace{1}(zCrossLoc(which));
              case {103,71} % 'g' key pressed
                D.isGood = 1;
                break;
              case {98,66} % 'b' key pressed
                D.isGood = 0;
                break
              case 113 % 'q' key pressed
                       % quit with empty output
                %keyboard
                Summary = struct;
                MEPs = struct;
                close(fighandle);
                return;
              case {44,60}
                % ,/< key pressed
                % go back one trial
                if i>1
                    skipSwitch = -1;
                    break;
                end
                
              case {46,62}
                % ./> key pressed
                % go forward one trial, but only if the current
                % trial is scored
                try
                    assert(ismember(Summary.isGood(i),[0 1]));
                    skipFlag = 1;
                    break;
                catch
                    Summary.isGood
                end
                
              otherwise
                % display what key was actually pressed; handy for debugging
                disp(button);
            end
            
            % Re-Plot my stuff:
            axes(wholeTraceAxes);
            cla
            plot(M.traceWhole{1});
            plot(D.peakTime+mepTime(1),D.peak,'R*');
            plot(D.troughTime+mepTime(1),D.trough,'G*');
            
%             if bgAxes
%                 axes(bgAxes)
%                 cla
%                 bar(1,D.pre_act_trial_avg,'b');
%                 errorbars(1,D.pre_act_trial_avg,D.pre_act_trial_std);
%                 bar(2,D.peak-D.trough,'r');
%             end    
    
            axes(mepAxes);
            cla
            plot(M.trace{1});
            plot(D.peakTime,D.peak,'R*');
            plot(D.troughTime,D.trough,'G*');
            
        end
    else, D.isGood = 1; end
    
    M.isGood = D.isGood;
    D.amplitude = D.peak-D.trough;
    
    if ~skipSwitch
        if i==1 || length(Summary.amplitude)<i
            Summary = addstruct(Summary,D);
            MEPs = addstruct(MEPs,M);
        else
            Summary = setrow(Summary,i,D);
            MEPs = setrow(MEPs,i,M);
        end    
    
        i=i+1;
    elseif skipSwitch==1
        i=i+1;
    else
        i=i-1;
    end
    
end

close(fighandle);
