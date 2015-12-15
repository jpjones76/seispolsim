% Test routine xcalign.m for accuracy

% Uses a delta function inserted into noise at random times.
ntr = 1000;
f = 0;
n = 0;
g = 0;

for k = 1:1:ntr
    X = randn(1024,6);
    for j=1:6; c=ceil(rand*1024); X(c,j)=1000; end
    
    % Call xcalign
    [Y,t,c] = xcalign(X,length(X));
    
    % Do the circular shifts match up?
    [i1,j1] = find(X==1000);
    m = round(mean(i1));
    d = -t - i1 + repmat(m,size(t));
    if max(d) > 1
        f = f+1;
        disp(['Trial ' num2str(k) ', shifts varied from mean ' ...
              'delta function position by up to' num2str(max(d)) ]);
        disp(['Trial ' num2str(k) ', mean CC was only ' ...
            sprintf('%04.2f',mean(c))]);
    elseif max(d) == 1
        n = n+1;
    else
        g = g+1;
    end
end
disp(['Exact success rate: ' sprintf('%5.1f', 100*g/ntr) '%']);
disp(['Exact and off-by-1: ' sprintf('%5.1f', 100*(g+n)/ntr) '%']);
disp(['Failure rate: ' sprintf('%5.1f', 100*f/ntr) '%']);