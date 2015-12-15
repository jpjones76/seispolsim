function setfigdefs(varargin)
% setfigdefs(OPTIONS);
% setfigdefs('only',OPTIONS);
% setfigdefs('check');
%
%   Set figure defaults. Specify OPTIONS as 'param',value.
%
% SPECIAL BEHAVIORS
% varargin{1} = 'only'   Only set things explicitly stated in OPTIONS
% varargin{1} = 'check'  List properties, don't actually set anything
% 
% ____________________________________________________________________
% OPTIONS
% Name    Default   Type Other Values          Property Set
% ----    -------   ---- ------------          ------------
% afont   Bookman   str  type 'listfonts'      DefaultAxesFontname
% afw     bold      str  light, demi, normal   DefaultAxesFontWeight
% axfs    12        sca  n > 0                 DefaultAxesFontSize
% box     on        str  off                   DefaultAxesBox
% fr      painters  str  zbuffer, OpenGL       DefaultFigureRenderer
% lw      1         sca  n > 0                 DefaultLineLineWidth
% ms      8         sca  n > 0                 DefaultLineMarkerSize
% np      replace   str  add, replacechildren  DefaultAxesNextPlot
% td      in        str  in, out               DefaultAxesTickDir
% ti      latex     str  tex, none             DefaultTextInterpreter
% tl      0.005     sca  n > 0                 DefaultAxesFontWeight
% txfs    12        sca  n > 0                 DefaultTextFontSize
% tfw     bold      str  light, demi, normal   DefaultTextFontWeight
% ____________________________________________________________________
% ====================================================================
% 
% WARNING: Functionality limited since Matlab 2014b broke setting the
% default interpreter (and everything else, really)
% 
% =======================================================================
% Author: Joshua Jones, highly.creative.pseudonym@gmail.com
% Version: 1.1, 2015-12-15


afont = 'Bookman';
afw = 'bold';
axfs = 12;
box = 'on';
check = 0;
lw = 1;
lms = 8;
fr = 'painters';
np = 'replace';
only = 0;
td = 'in';
tfont = afont;
tfw = afw;
ti = 'latex';
tl = 0.005;
txfs = 12;

if nargin
    j = 1;
    if strcmpi(varargin{1},'only')
        only = 1;
        j = j+1;
        strings = {};
    elseif strcmpi(varargin{1},'check')
        check = 1;
    end
    while j < nargin
        eval([varargin{j} '= varargin{j+1};']);
        if only
            strings = [strings; {varargin{j}}];
        end
        j = j+2;
    end
end

if only | check
    params = ...
        [{'box'}     {'DefaultAxesBox'}         {'box'};     ...
         {'axfs'}    {'DefaultAxesFontSize'}    {'axfs'};    ...
         {'afw'}     {'DefaultAxesFontWeight'}  {'afw'};     ...
         {'afont'}   {'DefaultAxesFontName'}    {'afont'};   ...
         {'np'}      {'DefaultAxesNextPlot'}    {'np'};      ...
         {'td'}      {'DefaultAxesTickDir'}     {'td'};      ...
         {'tl'}      {'DefaultAxesTickLength'}  {'[tl tl]'}; ...
         {'fr'}      {'DefaultFigureRenderer'}  {'fr'};      ...
         {'lw'}      {'DefaultLineLineWidth'}   {'lw'};      ...
         {'lms'}     {'DefaultLineMarkerSize'}  {'lms'};     ...
         {'tfont'}   {'DefaultTextFontName'}    {'tfont'};   ...
         {'txfs'}    {'DefaultTextFontSize'}    {'txfs'};    ...
         {'tfw'}     {'DefaultTextFontWeight'}  {'tfw'};     ...
         {'ti'}      {'DefaultTextInterpreter'} {'ti'}];
    
    if check
        for k = 1:1:size(params,1);
            tmp = get(0,params{k,2});
            if ~isstr(tmp)
                tmp = num2str(tmp);
            end
            disp([params{k,2} ' = ' tmp]);
        end
    else
        for j = 1:1:numel(strings)
            for k = 1:1:size(params,1);
                if strcmpi(params{k,1},strings{j})
                    eval(['set(0,params{k,2},' params{k,3} ');']);
                end
            end
        end
    end
else
    set(0, ...
        'DefaultAxesBox',         box,     ...
        'DefaultAxesFontSize',    axfs,    ...
        'DefaultAxesFontWeight',  afw,     ...
        'DefaultAxesFontName',    afont,   ...
        'DefaultAxesNextPlot',    np,      ...
        'DefaultAxesTickDir',     td,      ...
        'DefaultAxesTickLength',  [tl tl], ...
        'DefaultFigureRenderer',  fr,      ...
        'DefaultLineLineWidth',   lw,      ...
        'DefaultLineMarkerSize',  lms,     ...
        'DefaultTextFontName',    tfont,   ...    
        'DefaultTextFontSize',    txfs,    ...      
        'DefaultTextFontWeight',  tfw,     ...
        'DefaultTextInterpreter', ti );    
end

