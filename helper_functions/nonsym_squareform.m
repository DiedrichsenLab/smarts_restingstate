function [x] = nonsym_squareform(x,varargin)
% Reformats a matrix between a vector and a square form.
% [x] = nonsym_squareform(x,varargin)
%
% INPUT:
%       x:  1xN vector to convert into a matrix, or
%           NxN to convert into a vector. The size of x is determined by
%           dim
% OUTPUT:
%       x:      converted vector/matrix
% VARARGIN:
%   'dim':  'full' or 'red'
%           estimates the true size of the returned vector/matrix.
%           'full' keeps/assumes diagonal elements are present,
%           'red'  removes/assumes diag. elements are not present
% Naveed Ejaz (ejaz.naveed@gmail.com)
%
% Also see, squareform for symmetric conversion

% 0. Parse inputs
dim = 'full';   % full/red
vararginoptions(varargin,{'dim'});

% 1. Estimate size of conversion
[nRow,nCol] = size(x);
N           = ceil(sqrt(nRow*nCol));    % true size of the matrix
diag_idx    = 1:(N+1):(N^2);            % location of the diagonal elements

% 2. Perform conversation
switch(dim)
    case 'full'
        if nRow > 1
            x = reshape(x',1,numel(x));
        else
            x = reshape(x,N,N)';
        end;
    case 'red'
        if nRow > 1
            try
                x = reshape(x',1,numel(x));
                x(diag_idx) = [];
            catch
                x(diag_idx) = [];
            end;
        else
            for i=diag_idx
                x = [x(1:i-1) 0 x(i:end)];
            end;
            x = reshape(x,N,N)';
        end;
end;
