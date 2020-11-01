function [x , f] = obngifft(k,ft,n,rflag)
% Fast Fourier INVERSE Transform of the pair (k,ft) into (x,f).  The length of 
% x and f must be an even number, preferably a power of two.  
% x starts with zero.
% If rflag is set then only the real part of f is returned.
N         = size(k,n);
x         = x_of_k_ng(k,n);
dx        = diff(x);
Period    = N*dx(1);
inull     = N/2;
if n == 2
    f     = (N/Period)*ifft(circshift(ft,-[0 inull-1]),[],n);
elseif n == 1    
    f     = (N/Period)*ifft(circshift(ft,-[inull-1 0]),[],n);
end
if (nargin == 4)
    if (rflag == '1')
        f = real(f);
    end
end
