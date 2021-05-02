function [k , ft] = obngfft(x,f,n)
% Fast Fourier Transform of the pair (x,f) into (k,ft).  The length of 
% x and f must be an even number, preferably a power of two.  The index of
% the zero mode is inull=N/2. fft is performed along the nth dimension of
% f.
N         = size(x,n);
k         = k_of_x_ng(x,n);
dx        = diff(x,[],n);
Period    = N*dx(1);
inull     = N/2;
if n == 2
    ft    = (Period/N)*circshift(fft(f,[],n),[0 inull-1]);
elseif n == 1
    ft    = (Period/N)*circshift(fft(f,[],n),[inull-1 0]);
end

