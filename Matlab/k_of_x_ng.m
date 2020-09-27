function k = k_of_x_ng(x,n)
% Wave vector from x-vector
%
N      = size(x,n);
diffx  = diff(x,[],n);
dx     = diffx(1,1);
dk     = (2*pi)/(N*dx);
inull  = N/2;
k      = dk*(linspace(1,N,N)-inull);