%% compute the first nzeros zeros of the bessel functions of the
%% norders first orders of the nkind-th kind and saves the result in the
%% current directory in dht.mat. Needs the routine bessel_zeros.m by Adam
%% Wyatt (freely available in Matlab Central). The variable called 
%% precision is the convergence criterion for the algorithm.

clear all
tic
norders = 1024;
nzeros = 1024;
nkind = 1;
precision = 1e-15;

% scaling of the radial modes
c = zeros(norders/2+1,nzeros+1);
c(1,:) = bessel_zeros(nkind,0,nzeros+1,precision);
c(2:norders/2+1,1);
for n = 1:norders/2
    c(n+1,2:nzeros+1) = bessel_zeros(nkind,n,nzeros,precision);
end
toc
save('dht.mat','c')

beep