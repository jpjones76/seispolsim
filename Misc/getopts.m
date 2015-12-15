function opts = getopts(varargin)
% opts = getopts(varargin)
%   Parse input arguments in 'name', value pairs, return values in
% structure opts as opts.name  =  value. Also accepts structures as input
% arguments; each variable in each input structure is copied into the
% variable of the same name in opts. 
%   In general, it's a good idea to pass existing options to getopts with
% opts = getopts(opts, 'name', value)
%
% EXAMPLE:
% opts = getopts('fs', fs, 'L', L);
% 
% ======================================================================
% Author: Joshua Jones, highly.creative.pseudonym@gmail.com
% Version: 1.1, 2015-10-12

j = 1;
if isstruct(varargin{1})
    opts = varargin{1};
    j = 2;
end
while j <= nargin
    if isstruct(varargin{j})
        S = varargin{j};
        F = fieldnames(S);
        for f = 1:1:numel(F)
            opts.(F{f}) = S.(F{f});
        end
        j = j+1;
    else
        opts.(varargin{j}) = varargin{j+1};
        j = j+2;
    end
end
