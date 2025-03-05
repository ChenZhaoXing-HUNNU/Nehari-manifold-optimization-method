function y = dst2(x, mrows, ncols)
%DST2 Two-dimensional discrete sine Transform.
%   DST2(X) returns the two-dimensional sine transform of matrix X.
%   If X is a vector, the result will have the same orientation.
%
%   DST2(X,MROWS,NCOLS) reshape matrix X to size MROWS-by-NCOLS before
%   transforming.
%
%   See also FFT, IFFT, DST, IDST, IDST2.

%   Copyright: Wei Liu, 2016.


if nargin==1
    if min(size(x)) == 1
        y = dst(x);
    else
        y = dst(dst(x).').';
    end
else
    y = reshape(dst2(reshape(x,[mrows,ncols])),size(x));
end
