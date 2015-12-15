function updatepr(t,T)
% updatepr(n,N)
% 
% (NO OUTPUTS)
%   Updates progress on STDOUT as a percentage at iteration #n of N.
% =======================================================================
% Author: Joshua Jones, highly.creative.pseudonym@gmail.com
% Version: 0, written before 2010-01-01

if t == 0
    fprintf(1, ...
            [datestr(now,31) ': Starting ' num2str(T) ...
            ' calculations...\n']);
    fprintf(1,'Progress: %05.02f%c',0,char(37));
elseif t == T
    fprintf(1,'\b\b\b\b\b\b%05.02f%c. Done!\n',100.0,char(37));
else
    fprintf(1,'\b\b\b\b\b\b%05.02f%c',100*t/T,char(37));
end
