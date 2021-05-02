function x = x_of_k_ng(k,n)
% x-vector with leading zero from k-vector 
%
N      = size(k,n);
diffk  = diff(k,1,n);
dk     = diffk(1,1);
dx     = 2*pi/(N*dk);
x      = dx*(linspace(1,N,N)-1);