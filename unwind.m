function u = unwind(z)
%UNWIND  Unwinding number.
%   UNWIND(A) is the (scalar) unwinding number.

%   Reference: M. Aprahamian and N. J. Higham.
%   The matrix unwinding function, with an application to computing the
%   matrix exponential.  SIAM J. Matrix Anal. Appl., 35(1):88-109, 2014.

%   Mary Aprahamian and Nicholas J. Higham, 2013.

u = ceil( (imag(z) - pi)/(2*pi) );
