function x = idst2(y,mrows,ncols)
%IDST2 Two-dimensional inverse discrete sine transform.
%   IDST2(Y) returns the two-dimensional inverse sine transform of matrix
%   Y.  If Y is a vector, the result will have the same orientation.
%
%   IDST2(Y,MROWS,NCOLS) reshape matrix Y to size MROWS-by-NCOLS before
%   transforming.
%
%   See also FFT, IFFT, DST, IDST, DST2.

%   Copyright: Wei Liu, 2016.


if nargin == 1
    [mrows,ncols] = size(y);
end

x = 2^2/(mrows+1)/(ncols+1)*dst2(y,mrows,ncols);
