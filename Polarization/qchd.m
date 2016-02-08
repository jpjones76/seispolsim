function d = qchd(h1, h2, g)
% d = qchd(h1, h2, g);
%
% Quadratic-chi-square histogram distance (QCHD with m = 0.5) [1]; returns
% distances between matching column indices in h1, h2. 
%
% Reference: Pele, O., & Werman, M. (2010). The quadratic-chi histogram
% distance family. In Computer Vision - EECV 2010 (pp. 749-762). Springer
% Berlin Heidelberg. 
% 
% ========================================================================
% Author: Joshua Jones, highly.creative.pseudonym@gmail.com
% Version: 1.0, 2016-01-27

Z = (h1+h2)'*g;
Z(Z==0) = 1;
Z = (h1-h2)'./ sqrt(Z);
d = min( 1 , max( 0 , sqrt( max( 0 , sum( (Z*g).*Z , 2 ) ) ) ) );
