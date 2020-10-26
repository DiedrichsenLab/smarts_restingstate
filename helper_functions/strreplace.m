function stringout = strreplace(stringin,replaceme,replacewith)

% function stringout = strreplace(stringin,replaceme,replacewith)
%
% Similar to python's str.replace() function

%
%  Takes a string as input and an optional character that forms the
%  boundary between different "words", removes any instance of
%  the boundary and returns the remainder of the string as a cell
%  array of smaller strings.  If no characters are given, strsplit
%  uses a space as default.
%
%  EG: strsplit('Mary had a little lamb') would return:
%   { 'Mary' 'had' 'a' 'little' 'lamb' }

% JS -- 9/29/2006

%if nargin < 2
%  splitchar = ' ';
%end

if iscell(stringin)
  for string = 1:length(stringin(:))
    stringout{string} = strreplace(stringin{string},replaceme, ...
                                   replacewith);
  end
else
  
SplitByThese = strfind(stringin,replaceme);

if isempty(SplitByThese)
  stringout = stringin;
  return
end % else is assumed by return statement following if

cellsout = {};
%keyboard
while ~isempty(SplitByThese)
  cellsout{end+1} = stringin(1:SplitByThese(1)-1);
  stringin = stringin((SplitByThese(1)+length(replaceme)):end);
  SplitByThese = SplitByThese(2:end) - (SplitByThese(1)+length(replaceme)-1);
end

if ~isempty(stringin)
  cellsout{end+1} = stringin;
else
  cellsout{end+1} = '';
end


stringout = [];
for i = 1:length(cellsout)-1
  
  stringout = [stringout cellsout{i} replacewith];
  
end
stringout = [stringout cellsout{end}];

end