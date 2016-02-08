function G = gengdm(Nb,varargin)
t1 = floor(Nb/4);
t2 = floor(Nb/6);
c = [1 0 0 0 0];
wstr = {'hann';'gauss';'hann';'gauss';'gauss'};

if numel(varargin) > 0
    if ~isempty(varargin{1})
        t1 = varargin{1};
    end
    if numel(varargin) > 1
        if ~isempty(varargin{2})
            t2 = varargin{2};
        end
        if numel(varargin) > 2
            if ~isempty(varargin{3})
                c = varargin{3};
            end
            if numel(varargin) > 3
                if ~isempty(varargin{4})
                    wstr = varargin{4};
                end
            end
        end
    end
end

G = gdm(Nb, [t1 t2 t1 t2 t2], c, wstr);