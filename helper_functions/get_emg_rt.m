function [RT,avg_rt] = get_emg_rt(data_file,channels)
% function to load EMG data recorded from Signal 4.0 for visual inspection
% and manually marking EMG RT data
%
% Input: 
% data_file -- the full path of .mat file exported from Signal software
% channels -- an array of channel number to look at (can include a
% reference channel for artifect)
%
% Output:
% rt = array of EMG RT data for the entire data set, negative number
% indicates bad trials
%
% Usage:
% Move curser around to click the EMG-onset RT.  Click left-hand side of
% the figure to record a bad trial (NaN in the RT data array).

pretime = 500;
if ~exist(data_file,'file')
    [filename,path] = uigetfile('*.mat*','Please select a file');
    data_file = fullfile(path,filename);
end
data = load(data_file);

f = fields(data);
data = getfield(data,f{1});

figure('Position',[100,100,1200,800]);
% loop through total number of trials
% plotting all channels, but assuming only one channel to check for RT
RT = zeros(data.frames,1).* NaN;
n_channels = length(channels);
% n_channels = double(data.chans);
for t = 1:data.frames
    for ch = 1:n_channels
        subplot(n_channels,1,ch); box on; hold on;
        plot(data.values(:,channels(ch),t));
        ylim([-0.2 0.2]);
        line([1 data.points],[0.05 0.05],'Color','r');
        line([1 data.points],[0 0],'Color','k');
        line([1 data.points],[-0.05 -0.05],'Color','r');
        line([pretime pretime],[-0.2 0.2],'Color','k');
        ylabel('Volts');
        if ch == 1
            title(['Trial ' num2str(t)]);
        end
    end
    xlabel('Time (ms)');
    x = 0; y=0;
    [x y] = ginput(1);
    if x > 0
        RT(t) = round(x)-pretime;  % subtract pre_time
    end
    cla;
end
close gcf;
avg_rt = nanmean(RT(16:end));
save([data_file(1:end-4) '_EMGRT.mat'],'RT');

    