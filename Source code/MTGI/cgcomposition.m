
function [score,window]= cgcomposition(seq,varargin)
% Calculate CG content of the input data
% user's setting. 
%   seq-> Input genome


%   MTGIpick enables robust identification of genomic islands from a single genome
%   Qi Dai, 20 Apri 2016 



if isstruct(seq)
    seq = seqfromstruct(seq);
end

window = ceil(length(seq)/20);% ?cgthreshold = 0.5;
if nargin > 1
    if rem(nargin,2) == 0
        error('Bioinfo:IncorrectNumberOfArguments',...
            'Incorrect number of arguments to %s.',mfilename);
    end
    okargs = {'window','cgthreshold'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,numel(pname)));
        if isempty(k)
            error('Bioinfo:UnknownParameterName',...
                'Unknown parameter name: %s.',pname);
        elseif length(k)>1
            error('Bioinfo:AmbiguousParameterName',...
                'Ambiguous parameter name: %s.',pname);
        else
            switch(k)
                case 1   % window
                    window = pval;
                case 2  % gcthreshold
                    cgthreshold = pval;
            end
        end
    end
end

if window < 5 % window has a minimum length of 5
    window = 5;
elseif rem(window,2) == 0 % window must be odd
    window = window + 1;
end

if ischar(seq)
    seq = nt2int(seq);%
end

% Pad the end of the sequence to accomodate the window size
seq = [seq zeros(1,window-1)];
len = length(seq);

% Calculate nucleotide densities with a moving window
b = 1/window*ones(1,window);
a = 1;
sa = filter(b,a,seq==1);
sc = filter(b,a,seq==2);
sg = filter(b,a,seq==3);
st = filter(b,a,seq==4);

% Deal with the fact that filter starts and ends with zeros.
sa(1:window-1) = sa(1:window-1) * window./(1:window-1);
sc(1:window-1) = sc(1:window-1) * window./(1:window-1);
sg(1:window-1) = sg(1:window-1) * window./(1:window-1);
st(1:window-1) = st(1:window-1) * window./(1:window-1);

sa(len-window+2:len) = sa(len-window+2:len) * window./(window-1:-1:1);
sc(len-window+2:len) = sc(len-window+2:len) * window./(window-1:-1:1);
sg(len-window+2:len) = sg(len-window+2:len) * window./(window-1:-1:1);
st(len-window+2:len) = st(len-window+2:len) * window./(window-1:-1:1);

% Reduce the resulting vector to the correct subset
sa = sa(1+floor(window/2):len-floor(window/2));
sc = sc(1+floor(window/2):len-floor(window/2));
sg = sg(1+floor(window/2):len-floor(window/2));
st = st(1+floor(window/2):len-floor(window/2));
score=round((sc+sg)*window);
end
