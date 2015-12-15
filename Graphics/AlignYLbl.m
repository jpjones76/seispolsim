function AlignYLbl(varargin)
% AlignYLbl(fig,f,astr)
%
%   Align the Y labels of all axes in figure fig to fraction f of the
% normalized axes width. Set x alignment of Y labels to astr.
%   f and astr are optional.
% 
% DEFAULTS
% f = 0.07;
% astr = 'right';
% 
% =======================================================================
% Author: Joshua Jones, highly.creative.pseudonym@gmail.com
% Version: 1.1, last modified 2015-12-14

g = gcf;
f = .07;
hstr = 'right';
vstr = 'middle';
if nargin > 0
    g = varargin{1};
    if nargin > 1
        f = varargin{2};
        if nargin > 2
            hstr = varargin{3};
            if nargin > 3
                vstr = varargin{4};
            end
        end
    end
end
h = findall(g,'type','axes');
for i1 = 1:numel(h)
    if ~strcmpi(get(h(i1),'Tag'),'colorbar')
        h1 = get(h(i1),'ylabel');
        p1 = get(h1,'position');
        p2 = get(get(h1,'parent'),'xlim');
        p0 = -f*(p2(2)-p2(1));
        set(h1, ...
        'position',[p0+p2(1) p1(2) p1(3)], ...
            'verticalalignment', vstr,...
            'horizontalalignment', hstr);
    end
end
