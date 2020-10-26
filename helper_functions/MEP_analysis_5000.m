function [summ,meps] = MEP_analysis(varargin)
%varargin is sampling rate I input (either 1000 for wrong ones or 5000 for
%good ones)
% OK, first things first: process ALL .mat files in the directory,
% saving them as "filtered.mat" files.
% ---changed 07/2016 by Meret Branscheidt
% ---default options for the Palas project 2016

[FilePath] = uigetdir('Please select a directory');
cd(FilePath);

path = FilePath;
channel = 3;
gain = 200;
vararginoptions(varargin);
origDir = FilePath;

% cd(path)

summ = struct; meps = struct;
matSuff = '.mat';
filtSuff = '_filtered.mat';

allFiles = localFileSelector(matSuff);
filtFiles = localFileSelector(filtSuff);
origFiles = setdiff(allFiles,filtFiles);

% if ~isempty(filtFiles)g
%     if yesorno('Some files have been scored already.  Ignore those?');
%         useMe = ~ismember(strreplace(origFiles,matSuff,''),...
%                           strreplace(filtFiles,filtSuff,''));
%         origFiles = origFiles(useMe);
%     end
% end

for f=1:length(origFiles);
    %try % in case other .mat files exist (EMG?)
%         [S,M] = process_mep2(origFiles{f},...
        [S,M] = process_mep2(origFiles{f},...
                            struct('gain',gain,'channel',channel));
        save(strreplace(origFiles{f},matSuff,filtSuff),'S','M');
    %end
end

% 

cd(origDir);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local subroutines go below here

function filenames = localFileSelector(suffix,usenone)

if ~exist('usenone','var')
    usenone=0;
end

filenames = dir(['*' suffix]);
filenames = {filenames.name};

if usenone
    filenames(end+1) = {'<None>'};
end

return