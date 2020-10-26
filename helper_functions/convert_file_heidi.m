% 0. Setting up file name
fname       = 'rCST_aCST_lesioned.txt';          % name of file
fid         = fopen(fname,'r');                 % open file in read-only mode
outfname    = 'formatted_file.mat';             % name of formatted output file

% 1. Setting up header
D               = [];
Di              = [];
Di.SubjectID    = [];
Di.VisitWeek    = [];
Di.HemiTested   = [];
Di.MuscleTested = [];
Di.Outcome      = [];
Di.Trials       = [];

% 2. Main loop to read all lines in the file
lCount = 0;
while (true)
    line = fgets(fid);          % get the current line
    if line == -1               % reached end of the file
        break;
    end;
    lCount = lCount + 1;        % increment the counter
    
    % check if it is a line of interest (even numbered line)
    if lCount<=2 || mod(lCount,2)~=0
        continue;
    else
        % do data parsing into structure here
        disp(line);
        t = textscan(line,'%s%s%s%s%s%s','whitespace','|');     % read everything as string in first instance
        Di.SubjectID    = {strtrim(t{1}{1})};
        Di.VisitWeek    = str2double(strtrim(t{2}{1}));
        Di.HemiTested   = {strtrim(t{3}{1})};
        Di.MuscleTested = {strtrim(t{4}{1})};
        Di.Outcome      = {strtrim(t{5}{1})};
        Di.Trials       = str2double(strtrim(t{6}{1}));
        D               = addstruct(D,Di);
    end;
end;

% 4. Saving formatted data as file
save(outfname,'-struct','D');