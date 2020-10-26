function varargout = process_mep_premove(dd,flags)

% USAGE: [Summary,MEPs,flags] = process_mep_premove(dd,flags);
%
% INPUT:
% DataIn can be either a filename (full path), a datastructure
% (loaded from a file exported from Signal 4.0X), empty, or not
% supplied.  If empty or not supplied, user will be prompted to
% select a *.mat file (exported from Signal 4.0X).
%
% flags is a structure for analysis flags.  If not supplied, a
% quick dialog will appear to set the values.  If the structure is
% supplied, any missing fields will be filled with the default
% values, but no GUI will appear.  Fields are as follows:
%
%   channel -- data channel to analyze.  Default is -1, which
%      prompts the user to specify manually.
%   multiplier -- Value that needs to be multiplied to raw traces
%      to get values in mV.  Default is 1 (ie, values are in mV
%      already).
%   gain -- This is a legacy special case for the Viking EMG filter
%      in Pablo Celnik's lab at Johns Hopkins University, where the
%      output analog EMG values recorded by Signal are scaled by a
%      value called "gain" on the Viking console.  This is not used
%      by default, but if set then it will set the multiplier to
%      gain/200, overriding any value passed to multiplier.
%   background_noise -- boolean, to measure background noise.
%      Default is 1.
%   mepTime -- vector of times to use for recording MEPs.  Times
%      are specified as indices into the channel matrix.  Default
%      is [800:1100], which works for the standard Celniklab TMS
%      configurations but may not be appropriate for all configs.
%   preTime -- vector of times to use for measuring preactivation
%      of EMG.  Units as mepTime, default is [51:200].  Again, this
%      default works for the standard Celniklab TMS configurations,
%      but may not be appropriate for all configs.
%   totalTime -- vector of times to use for saving the MEP trace.
%      Units as for previous two flags.  Default is [1:750].
%      Again, this default works for the standard Celniklab TMS
%      configurations, but may not be appropriate for all configs.
%   handfilter -- boolean, to allow user to filter through MEP
%      traces using a GUI.  Default is 1 (highly recommended!)
%   active_cursor -- Value from MVC measures, used for drawing lines for 
%      active MEP plots at 20% MVC.
%   frames -- a vector with the frame numbers to be marked.  This requires
%      pre-existing processed data.  If empty, go through all the frames.
%
% GUI commands (only if handfilter is true):
%  q key -- Quit, returning empty output
%  g key -- Mark the current MEP as "good"
%  b key -- Mark the current MEP as "bad"
%  h key -- Bring up a quick "Help" dialog.
%  < key -- Return to the previous frame (ie, to re-score)
%  > key -- Advance to the next frame (only if current frame is
%           scored as good/bad)
%  Left mouse button -- Adjust the peak (Red *) to the nearest 
%           local minima/maxima to the cursor
%  Right mouse button -- Adjust the trough (Green *) to the nearest
%           local minima/maxima to the cursor
%  Middle mouse button -- Bring up a dialog to mark the trial good
%           or bad
%
% OUTPUT:
% Summary will be a structure, with every field containing a column
%   of data.  This is also saved as a text file and can be opened in
%   a spreadsheet program if you so choose.
% MEPs will have some of the same data as Summary, as well as MEP
%   traces themselves.  This will also be a column-oriented
%   structure; the MEP traces will be cell arrays.
% flags as output will store a copy of the option flags used.
%
% If no output is specified, the user will be prompted to save the
% data, and can then specify a filename into which the Summary and
% MEPs structures will be saved.


%%%%%%%%%%%
% This is essentially a streamlined version of mep_any, by 
%   Joe Galea
%   -- John Schlerf, August 2010
% Comments and help text updated
%   -- John Schlerf, September 2011
% More help updating, some streamlining to make it less
%   CelnikLab-focused.
%   -- John Schlerf, October 2011
% Updated for the default setup for SMARTS TMS measures
%   -- Jing Xu, Feb. 2013
% Modefied for the premovement data processing
%   -- Jing Xu, Apr. 2013
%%%%%%%%%%%

% load both chennels at the bottom, and one channel on top for MEP
% processing.
% check the state (in title)
% channel 1: mark RT & MEP
%            - if inside RT, mark RT and NaN
% channel 2: mark MEP & mirror movement RT
% good trial:
%   - no EMG activitiy before TMS pulse
%   - MEP is before RT

% set up default option flags:
% parameters for JHU data
% delay = 0;
% parameters for Vitznau with 312ms delay due to digital EMG input
delay = 0;
optflags = struct('channel',[2 3],'multiplier',1,'background_noise',1,...
                  'mepTime',[500:1750]+delay,'preTime',[500:750]+delay,...
                  'totalTime',[1:5000],'gain',[],'frames',[]);

% parameters for UZ
% delay = 320;
% optflags = struct('channel',[2 3],'multiplier',1,'background_noise',1,...
%                   'mepTime',[1850:3000]+delay,'preTime',[1850:2000]+delay,...
%                   'totalTime',[1:5000],'gain',[],'frames',[]);


% optflags = struct('channel',[2 3],'multiplier',1,'background_noise',1,...
%                   'mepTime',[2170:3000],'preTime',[2170:2320],...
%                   'totalTime',[1:5000],'gain',[],'frames',[]);
% parameters for CUP1002
% optflags = struct('channel',[2 3],'multiplier',1,'background_noise',1,...
%                   'mepTime',[1850:2700],'preTime',[1850:2000],...
%                   'totalTime',[1:5000],'gain',[],'frames',[]);

% parse my input:
if exist('dd','var')
    if ~exist('flags','var') 
        if isstruct(dd) && any(ismember(fields(dd),fields(optflags)))
            flags = dd;
            dd = [];
        end
    end
end

if ~exist('dd','var') || isempty(dd)
    disp('Load data');
    [file] = find_file('IHIR*.mat');
    if isempty(file)
        disp('No file is selected.');
        return;
    end
elseif ischar(dd)
    if exist(dd,'file')
        file = dd;
    else ischar(dd)
        disp('Cannot find specified file, please reselect.');
        [file] = find_file('*.mat*');
    end
end

if exist('flags','var')
    f = fields(flags);
    for i=1:length(f)
        optflags = setfield(optflags,f{i},getfield(flags,f{i}));
    end
else
    optflags = optflagsgui(optflags);
    if isempty(optflags)
        return;
    end
end

% load if needed:
if exist('file','var')
    dd = load(file); 
    [pathstr,name,ext]= fileparts(file);
    disp(['Processing: ' name]);
    f = fields(dd);
    dd = getfield(dd,f{1});
    mep_file = [file(1:end-4) '.MEP.mat'];
    if exist(mep_file,'file')
        load(mep_file);
    end
end

if isempty(optflags.frames) % take all the frames
    numFrames = dd.frames; % how many MEPs (Frames) do we have?
else
    numFrames = length(optflags.frames);
    if ~exist(mep_file,'file')
        error('err:mep_file','Can not load pre-processed MEP file.');
    end
end
% channels=dd.chans;
% good_trial=ones(1,numFrames); % we'll use this for filtering
% thresh = 0.25; % I don't remember how this is used

% Select the channel that contains the data:
data_channel = optflags.channel;
if data_channel == -1 % prompt for data_channel if not provided
    Dchannel = {''};
    while isempty(Dchannel{1})
        Dchannel = inputdlg('Please input the data channel.');
    end
    data_channel = eval(Dchannel{1});
end

% And the rest of the options, which don't need prompting:
multiplier = optflags.multiplier;
background_noise = optflags.background_noise;
mepTime = optflags.mepTime;
totalTime = optflags.totalTime;
preTime = optflags.preTime;
MEP_Total_Offset = mepTime(1)-totalTime(1);
% OK, legacy special case here:
if ~isempty(optflags.gain)
    multiplier = optflags.gain/200;
end

%% now process all data that I want, MEP by MEP

% Start by setting up a figure for the GUI, if we're using it:
scrsz = get(0, 'ScreenSize');
fighandle = figure('Position',scrsz); 
mepAxes = axes('position',[0.05,0.5,0.9,0.5]); % left, bottom, width, height
hold on;
wholeTraceAxes(1) = axes('position',[0.05,0.26,0.7,0.2]);
wholeTraceAxes(2) = axes('position',[0.05,0.05,0.7,0.2]);
hold on;
if background_noise
    bgAxes(1) = axes('position',[0.75,0.26,0.2,0.2]);
    bgAxes(2) = axes('position',[0.75,0.05,0.2,0.2]);
    hold on;
else
    bgAxes(1:2) = 0;
end

% Initialize my structures:
if ~exist('Summary','var')
    Summary = struct;
end
if ~exist('MEPs','var')
    MEPs = struct;
end

i=1;
while i <= numFrames
    if isempty(optflags.frames)
        f = i;
    else
        f = optflags.frames(i);
    end
    
    skipSwitch = 0;
    skipFlag = 0;
    
    % structs to hold current frame data
    M = struct;
    D = struct;
    if ~exist(mep_file,'file')
        % Pass frameinfo values into output:
        for c = 1:length(data_channel)
            D.channel(c) = dd.frameinfo(f);
            M.channel(c) = dd.frameinfo(f);
        end   
        % initial rows which will go to MEP structure:
        for c=1:length(data_channel)
            M.channel(c).trace = {(multiplier)*squeeze(dd.values(mepTime,data_channel(c),f))};
            M.channel(c).traceWhole = {(multiplier)*squeeze(dd.values(totalTime,data_channel(c),f))};
            M.channel(c).time = {mepTime};
            %keyboard;
        end

        % Guess min and max of each trial
        for c=1:length(data_channel)
            [D.channel(c).peak,D.channel(c).peakTime] = max(M.channel(c).trace{1});
            D.channel(c).peakTime = D.channel(c).peakTime+mepTime(1)-preTime(2); % start from time of TS pulse
            [D.channel(c).trough,D.channel(c).troughTime] = min(M.channel(c).trace{1});
            D.channel(c).troughTime = D.channel(c).troughTime+mepTime(1)-preTime(2);
            if c==2 && dd.frameinfo(f).state==1
                D.channel(c).peak = NaN;
                D.channel(c).trough = NaN;
                D.channel(c).peakTime = NaN;
                D.channel(c).troughTime = NaN;
            end
            D.channel(c).RT = NaN;
        end    

        % if we're collecting bg noise, do it here:
        if background_noise
            for c=1:length(data_channel)
                D.channel(c).rms = sqrt(mean(dd.values(preTime,data_channel(c),f).^2));
                D.channel(c).pre_act_trial_avg = mean(abs(dd.values(preTime, ...
                                                        data_channel(c),f)));
                D.channel(c).pre_act_trial_std = std(abs(dd.values(preTime, ...
                                                        data_channel(c),f)));
                M.channel(c).bg = {squeeze(dd.values(preTime,data_channel(c),f))};
                M.channel(c).bg_rec = {abs(M.channel(c).bg{1})};
            end
        end
    else  %load pre-processed data
        for c=1:length(data_channel)
            % Summary goes to D, MEPs goes to M
            fnames = fieldnames(Summary.channel);
            for k=1:length(fnames)
                D.channel(c).(fnames{k}) = Summary.channel(c).(fnames{k})(f);
            end
            fnames = fieldnames(MEPs.channel);
            for k=1:length(fnames)
                M.channel(c).(fnames{k}) = MEPs.channel(c).(fnames{k})(f);
            end
        end
    end
    
    % start plotting
    % whole traces
    for c=1:length(data_channel)
        axes(wholeTraceAxes(c)); cla; hold on;
        plot(M.channel(c).traceWhole{1});
        plot(D.channel(c).peakTime+MEP_Total_Offset,D.channel(c).peak,'R*');
        plot(D.channel(c).troughTime+MEP_Total_Offset,D.channel(c).trough,'G*');
        if c==1
            t=title('Entire TS trace (press H for Help)');
        else
            t=title('Entire CS trace');
        end
        % move title down a little:
        if max(abs(ylim))<1.0
            ylim([-1.0 1.0]);
        else
            ylim('auto');
        end
        YL = ylim(wholeTraceAxes(c));
        TL = get(t,'position');
        TL(2) = YL(2)-0.1*diff(YL);
        set(t,'position',TL);
        ylabel('mV');
        if c==2
            xlabel('Time (msec)');
        end
        % lines for RMS threshold
        XL = xlim(wholeTraceAxes(c));
        line(XL,[0 0],'Color','k');
        line(XL,[0.05 0.05],'Color','k');
        line(XL,[-0.05 -0.05],'Color','k');
        line([2000 2000]+delay,YL,'Color','k');
        if ~isnan(D.channel(c).RT)
            % draw a line for RT
            line([D.channel(c).RT D.channel(c).RT]+MEP_Total_Offset,ylim,'Color','k');
        end
    end
    % background bars
    if bgAxes(c)
        for c=1:length(data_channel)
            axes(bgAxes(c)); cla; hold on;
            bar(1,D.channel(c).pre_act_trial_avg,'b');
            errorbars(1,D.channel(c).pre_act_trial_avg,D.channel(c).pre_act_trial_std);
            bar(3,D.channel(c).peak-D.channel(c).trough,'r');
            ylabel('mV');
            set(gca,'xtick',[1,3]);
            set(gca,'xticklabel',{'BG','MEP'});
        end
    end

    % start marking
    for c=1:length(data_channel)
        axes(mepAxes); cla; hold on;
        % plot MEP, and adjust max/min if necessary:
        plot(M.channel(c).trace{1});
        plot(D.channel(c).peakTime,D.channel(c).peak,'R*');
        plot(D.channel(c).troughTime,D.channel(c).trough,'G*');
        if ~isfield(D.channel(c),'isGood') || isempty(D.channel(c).isGood)
            if c==1
                t=title(sprintf('TS MEP %d of %d; State: %d',f,dd.frames,dd.frameinfo(f).state));
            else
                t=title(sprintf('CS MEP %d of %d; State: %d',f,dd.frames,dd.frameinfo(f).state));
            end
        else
            if D.channel(c).isGood
                if c==1
                    t=title(sprintf('TS MEP %d of %d; State: %d;    Good',f,dd.frames,dd.frameinfo(f).state));
                else
                    t=title(sprintf('CS MEP %d of %d; State: %d;    Good',f,dd.frames,dd.frameinfo(f).state));
                end
            else
                if c==1
                    t=title(sprintf('TS MEP %d of %d; State: %d;    Bad',f,dd.frames,dd.frameinfo(f).state));
                else
                    t=title(sprintf('CS MEP %d of %d; State: %d;    Bad',f,dd.frames,dd.frameinfo(f).state));
                end
            end
        end
        % move axis title down a little:
        YL = ylim;
        TL = get(t,'position');
        TL(2) = YL(2)-0.1*diff(YL);
        set(t,'position',TL);
        % lines for RMS threshold
        XL = xlim(mepAxes);
        line(XL,[0 0],'Color','k');
        line(XL,[0.05 0.05],'Color','k');
        line(XL,[-0.05 -0.05],'Color','k');
        line([150 150],ylim,'Color','k');
        if ~isnan(D.channel(c).RT)
            % draw a line for RT
            line([D.channel(c).RT D.channel(c).RT],ylim,'Color','k');
        end

        while 1
            % Let the user click/press a key
            [x,y,button]=ginput(1);
            if x < 0
                x = NaN;
            else
                x=round(x);
            end
            if ~isempty(button) % this happens with F1, for example
              switch button % deal with clicks
                case {1} % left mouse button, MEP peak
                  % Use the nearest local max or min:
                  dt = [0; sign(diff(M.channel(c).trace{1}))];
                  zCrossLoc=find(diff(dt));
                  if ~isnan(x)
                      distSq = (zCrossLoc - x).^2;
                      [minDist,which] = min(distSq);
                      D.channel(c).peakTime = zCrossLoc(which);
                      D.channel(c).peak = M.channel(c).trace{1}(zCrossLoc(which));
                  else
                      D.channel(c).peakTime = NaN;
                      D.channel(c).peak = NaN;
                  end
                  D.channel(c).RT = NaN; % reset RT to NaN to erase previous markings
                case {2} % middle mouse button, EMG RT
                  D.channel(c).RT = x;    
                case {3} % right mouse button, MEP trough
                  % Use the nearest local max or min:
                  dt = [0; sign(diff(M.channel(c).trace{1}))];
                  zCrossLoc=find(diff(dt));
                  if ~isnan(x)
                      distSq = (zCrossLoc - x).^2;
                      [minDist,which] = min(distSq);
                      D.channel(c).troughTime = zCrossLoc(which);
                      D.channel(c).trough=M.channel(c).trace{1}(zCrossLoc(which));
                  else
                      D.channel(c).troughTime = NaN;
                      D.channel(c).trough = NaN;
                  end
                case 103 % 'g' key pressed, mark good:
                  D.channel(c).isGood = 1;
                  break;
                case 98 % 'b' key pressed, mark bad:
                  D.channel(c).isGood = 0;
                  break;
                case 113 % 'q' key pressed, quit with empty output:
                  Summary = struct;
                  MEPs = struct;
                  close(fighandle);
                  return;
                case 104 % 'h' key pressed; bring up "help" dialog:
                  uiwait(helpdlg({'List of GUI commands',...
                                  'MOUSE:',...
                                  '  Right Click -- Move Trough',...
                                  '  Left Click -- Move Peak',...
                                  '  Middle Click -- locate the corsor for RT',...
                                  '',...
                                  'KEYBOARD:'...
                                  '  G -- Mark MEP as Good',...
                                  '  B -- Mark MEP as Bad',...
                                  '  < -- Move backward one frame',...
                                  '  > -- Move forward one frame',...
                                  '  H -- Bring up this HELP menu',...
                                  '  Q -- Quit, returning nothing'},...
                                 'GUI Help'));
                case {44,60}
                  % ,/< key pressed
                  % go back one trial
                  if f>1
                      skipSwitch = -1;
                      skipFlag = 1;
                      break;
                  end
                case {46,62}
                  % ./> key pressed
                  % go forward one trial; but only if the current trial is scored
                  if isfield(Summary,'channel')
                      try
                          assert(ismember(Summary.channel(c).isGood(f),[0 1]));
                          skipSwitch = 1;
                          skipFlag = 1;
                          break;
                      catch
                          Summary.channel(c).isGood;
                      end
                  end
                otherwise
                  % display what key was actually pressed; handy for debugging
                  disp(button);
              end
            end
            
            % Re-Plot stuff since I've rescored:
            axes(wholeTraceAxes(c));
            cla; hold on;
            plot(M.channel(c).traceWhole{1});
            plot(D.channel(c).peakTime+MEP_Total_Offset,D.channel(c).peak,'R*');
            plot(D.channel(c).troughTime+MEP_Total_Offset,D.channel(c).trough,'G*');
            line(xlim,[0.05 0.05],'Color','k');
            line(xlim,[-0.05 -0.05],'Color','k');
            if ~isnan(D.channel(c).RT)
                % draw a line for RT
                line([x x]+MEP_Total_Offset,ylim,'Color','k');
            end
            
            if bgAxes(c)
                axes(bgAxes(c));
                cla; hold on;
                bar(1,D.channel(c).pre_act_trial_avg,'b');
                errorbars(1,D.channel(c).pre_act_trial_avg,D.channel(c).pre_act_trial_std);
                bar(2,D.channel(c).peak-D.channel(c).trough,'r');
            end    

            axes(mepAxes);
            cla; hold on;
            plot(M.channel(c).trace{1});
            plot(D.channel(c).peakTime,D.channel(c).peak,'R*');
            plot(D.channel(c).troughTime,D.channel(c).trough,'G*');
            line(xlim,[0.05 0.05],'Color','k');
            line(xlim,[-0.05 -0.05],'Color','k');
            line([150 150],ylim,'Color','k');
            if ~isnan(D.channel(c).RT)
                line([x x],ylim,'Color','k');
            end
            if ~isfield(D.channel(c),'isGood') || isempty(D.channel(c).isGood)
                if c==1
                    t=title(sprintf('TS MEP %d of %d; State: %d',f,dd.frames,dd.frameinfo(f).state));
                else
                    t=title(sprintf('CS MEP %d of %d; State: %d',f,dd.frames,dd.frameinfo(f).state));
                end
            else
                if D.channel(c).isGood
                    if c==1
                        t=title(sprintf('TS MEP %d of %d; State: %d;    Good',f,dd.frames,dd.frameinfo(f).state));
                    else
                        t=title(sprintf('CS MEP %d of %d; State: %d;    Good',f,dd.frames,dd.frameinfo(f).state));
                    end
                else
                    if c==1
                        t=title(sprintf('TS MEP %d of %d; State: %d;    Bad',f,dd.frames,dd.frameinfo(f).state));
                    else
                        t=title(sprintf('CS MEP %d of %d; State: %d;    Bad',f,dd.frames,dd.frameinfo(f).state));
                    end
                end
            end
            YL = ylim;
            TL = get(t,'position');
            TL(2) = YL(2)-0.1*diff(YL);
            set(t,'position',TL);
        end
        if skipFlag
            break;
        end
        M.channel(c).isGood = D.channel(c).isGood;
        D.channel(c).amplitude = D.channel(c).peak-D.channel(c).trough;
    end
    if ~skipSwitch
        if ~exist(mep_file,'file')
            if f==1
                Summary = addstruct(Summary,D);
                MEPs = addstruct(MEPs,M);
            else
                for c=1:length(data_channel)
                    if length(Summary.channel(c).isGood)<f  % adding a new row
                        Summary.channel(c) = addstruct(Summary.channel(c),D.channel(c));
                        MEPs.channel(c) = addstruct(MEPs.channel(c),M.channel(c));
                    else
                        Summary.channel(c) = setrow(Summary.channel(c),f,D.channel(c));
                        MEPs.channel(c) = setrow(MEPs.channel(c),f,M.channel(c));
                    end    
                end
            end
        else
            for c=1:length(data_channel)
                Summary.channel(c) = setrow(Summary.channel(c),f,D.channel(c));
                MEPs.channel(c) = setrow(MEPs.channel(c),f,M.channel(c));
            end
        end
        i=i+1;
    elseif skipSwitch==1
        i=i+1;
    else
        i=i-1;
    end
%    keyboard;
end

close(fighandle);

if nargout == 0
    % Prompt user to save data.
    yesno = questdlg({'Do you want to save your work?', ...
                      'Answering NO means that the data',...
                      'variables will be added to the',...
                      'base workspace.'},...
                     'Save','Yes','No','Yes');
    if isequal(yesno,'Yes')
        outfilename = file(1:end-4);
        outmatfile = [outfilename '.MEP.mat'];
        outtxtfile1 = [outfilename '.MEP_ch1.txt'];
        outtxtfile2 = [outfilename '.MEP_ch2.txt'];
        
%        [outmatfile,outmatpath] = uiputfile('*.mat','Save data as mat file:');
%        outmatfile = fullfile(outmatpath,outmatfile);
        save(outmatfile,'Summary','MEPs');
%        [outtxtfile,outtxtpath] = uiputfile('*.txt','Save summary as txt file:');
%        outtxtfile = fullfile(outtxtpath,outtxtfile);
        dsave(outtxtfile1,Summary.channel(1));
        dsave(outtxtfile2,Summary.channel(2));

%         uiwait(msgbox(sprintf(['Created a file called %s, ' ...
%                             'in folder %s.  The file ' ...
%                             'contains two variables: Summary ',...
%                             'and MEPs.'],outmatfile,outmatpath),'Status'));
%         uiwait(msgbox(sprintf(['Created a file called %s, ' ...
%                             'in folder %s.  The file ' ...
%                             'contains Summary.'],outtxtfile,outtxtpath),'Status'));
    else
        assignin('caller','Summary',Summary);
        assignin('caller','MEPs',MEPs);
    end
    
elseif nargout == 1
    varargout = {Summary};
elseif nargout == 2
    varargout = {Summary, MEPs};
else
    varargout = {Summary,MEPs,optflags};
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OF = optflagsgui(optflags)

F = fields(optflags);
for f=1:length(F)
    default{f} = prettystr(optflags.(F{f}));
end

newvals = inputdlg(F,'Analysis Options',1,default);

if isempty(newvals)
    OF = [];
else
    for f=1:length(F)
        if ~isempty(newvals{f})
            OF.(F{f}) = eval(newvals{f});
        else
            OF.(F{f}) = [];
        end
    end
end

return
