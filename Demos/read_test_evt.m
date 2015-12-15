function varargout = read_test_evt(evt, fs, fc, fl)
switch lower(evt)
    case 'tr1'
        % Tremor-like noise seen 2013 01/01, late at night
        fname = 'Hoadley20130101044634.csv';
    
    case 'fsj'
        % Mn = 4.2 near Ft. St. John, BC
        fname = 'LocalEQ/Hoadley20130527223643.csv';
        tP = 42.62;
        phi = -46.345;
        fh = 12;
        tos = 0;
        titlestr = '$$\mathbf{M_w}$$\bf{ = 4.2, 2013 05/28, }$$\mathbf{\Delta = 5.0^{\circ}}$$';
        
    case 'okh'
        % Sea of Okhhostk M = 8.3 2013 05/24 (UTC)
        fname = 'Teleseism/Hoadley20130523235232.csv';
        tP = 22.8;
        phi = -48.340;
        fc = 0.2;
        fh = 4.5;
        tos = 461;
        titlestr = '$$\mathbf{M_w}$$\bf{ = 8.3, 2013 05/23, }$$\mathbf{\Delta = 50.6^\circ}$$';
           
    otherwise
        % Mn = 3.1, 67 km SE of Banff
        fname = 'LocalEQ/Hoadley20130629061956.csv'; 
        tP = 31.6;
        phi = -169.8294;
        fh = 20;
        tos = 5;
        titlestr = '$$\mathbf{M_N}$$\bf{ = 3.1, 2013 06/29, }$$\mathbf{\Delta = 244 km}$$';
end
disp(['Reading ' fname '...']);
X = read_hoadley(['/Hoadley/Samples/' fname], fs, fc, fl, fh);
disp('Data loaded.');

% Remove stations 3, 11
X = X(1+10*fs:130*fs,[1:6 10:30 34:36]);
[~,Nc] = size(X); Nk = Nc/3;

% Align X in time based on Z component.
[~,sh] = xcalign(X(:,1:3:Nc));
for k = 1:1:Nk
    X(:,3*k-2:3*k) = circshift(X(:,3*k-2:3*k),[sh(k) 0]);
end

varargout{1} = X;
varargout{2} = tos;
if nargout > 2
    tP = varargout{3};
    phi = varargout{4};
    titlestr = varargout{5};
end