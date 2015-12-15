% Script to test seispol.m for accuracy
clear all; close all; fclose('all');

% These can be user-controlled
SN = 5;
av = 1;
dp = 0.01;
fca = 0;

%% Don't modify anything below here
% X spacing for plots
ex = .5*dp:dp:1-.5*dp;
ax = (ex*180*(1+fca))-90;

%% Generate synthetic (requires wavelet toolbox)
x1 = resample(wavefun('db32')', 1024, 16129);
x1 = [zeros(512,1); x1; zeros(512,1)];
t1 = find(abs(hilbert(x1))>0.1, 1 );
t2 = find(abs(hilbert(x1))>0.1, 1, 'last' );
x1 = x1./norm(x1);
px = sum(x1.^2);

%% First check: Is Az nearly due N? Is ellipticity low?
X = randn(2048,3);
X = X./repmat(sqrt(sum(X.^2,1)),2048,1);
pn = sum(X(:,2).^2);

% Control S:N 
x1 = x1.*(pn/px)*(10^(SN/20));
disp(['S/N = ' sprintf('%0.1f',10*log10(sum(x1.^2)/sum(X(:,2).^2))) ' dB']);
X(:,2) = X(:,2) + x1;
c = 2*max(max(abs(X)));
X = X./c;

[P,~,H] = seispol(X,'av',av,'h1',1,'dp',dp,'fca',fca);
testpolfig;
set(gcf,'Name','Az 0, In 90');

%% Next test: Az should be 45 degrees east of north
X = randn(2048,3);
X = X./repmat(sqrt(sum(X.^2,1)),2048,1);
X(:,2) = X(:,2) + x1;
X(:,3) = X(:,3) + x1;
c = 2*max(max(abs(X)));
X = X./c;
[P,~,H] = seispol(X,'av',av,'h1',1,'dp',dp,'fca',fca);
testpolfig;
set(gcf,'Name','Az 45, In 90');

%% Test 3: Az 45, In 60
X = randn(2048,3);
X = X./repmat(sqrt(sum(X.^2,1)),2048,1) + repmat(x1,1,3);
c = 2*max(max(abs(X)));
X = X./c;
[P,~,H] = seispol(X,'av',av,'h1',1,'dp',dp,'fca',fca);
testpolfig;
set(gcf,'Name','Az 45, In 60')