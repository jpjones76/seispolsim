% Script to check whether wavelet transform routines work correctly.
clear all;
close all;
fclose('all');

%% These can be altered
X = detrend(rand(2^16, 3), 'constant');
tol = 1.0e-9;
J = 4;
wlt = 'LA{16}';

%% Static parameters
E0 = sum(X(:).^2);
Er = {}; 
c = 1;

%% Check MODWPT
W = modwt(X, wlt, J, 1);
for j = 1:1:J
    En = zeros(size(X));
    for n = 1:2^j
        En = En + sum(W{j,n}(:).^2);
    end
    if abs(En - E0) > tol
        Er{c} = ['MODWPT is not preserving energy in scale ' num2str(j)];
        c = c+1;
    end
end
if isempty(Er)
    disp(['Energy preserved by each MODWPT scale j=1:' num2str(J) '.']);
end

%% Check inverse MODWPT
D = cell(size(W));
for j = 1:1:J; 
    Xr = zeros(size(X));
    for n = 1:1:2^j 
        D{j,n} = iwpt(W{j,n}, [j n], wlt, 1);
        Xr = Xr + D{j,n};
    end; 
    rerr = max(abs((X(:)-Xr(:))));
    if rerr > tol
        Er{c} = ['IMODWPT is not reconstructing in scale ' num2str(j) ...
            ' (err = ' num2str(rerr) ')' ];
        c = c+1;
    end
end
% Detail coeffs. don't preserve energy (Schwartz inequality)

if isempty(Er)
    disp(['Perfect reconstruction from each scale j=1:' num2str(J) '.']);
    disp('No errors found.');
else
    disp([num2str(c-1) ' errors found:']);
    for c = 1:1:numel(Er)
        disp(Er{c});
    end
end