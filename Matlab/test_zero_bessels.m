clear all
jmodes = 100
nth = 100
R = 5
% scaling of the radial modes
zero_bess_t = zeros(jmodes+1,nth);
%zero_bess_u = zeros(jmodes+1,nth);
zero_bess_v = zeros(jmodes+1,nth);
tic
for n = 0:nth/2
    %zero_bess_t(:,n+nth/2) = JnRoots(n,jmodes+1);
    %zero_bess_u(:,n+nth/2) = besselzero(n,jmodes+1,1);
    zero_bess_v(:,n+nth/2) = bessel_zeros(1,n,jmodes+1,1e-14);
    % WARNING: the 0-order bessel function of the 1st kind is described
    % by n-index = nth/2, etc. Therefore, expect a lot of confusion.
end

for n = 1:nth/2-1
    %zero_bess_t(:,n) = zero_bess_t(:,nth-n);
    %zero_bess_u(:,n) = zero_bess_u(:,nth-n);
    zero_bess_v(:,n) = zero_bess_v(:,nth-n);
end
%K_t = zero_bess_t/R;
%K_u = zero_bess_u/R;
K_v = zero_bess_v/R;

%plottestt = zeros(jmodes+1,nth);
%plottestu = zeros(jmodes+1,nth);
plottestv = zeros(jmodes+1,nth);
for n = -nth/2+1:nth/2
   % plottestt(:,n+nth/2) = besselj(n,zero_bess_t(:,n+nth/2));
    %plottestu(:,n+nth/2) = besselj(n,zero_bess_u(:,n+nth/2));
    plottestv(:,n+nth/2) = besselj(n,zero_bess_v(:,n+nth/2));
end
toc
%figure(1)
%pcolor(log10(abs(real(plottestt)))), shading flat
%colorbar
%figure
%pcolor(log(real(plottestu))), shading flat
%colorbar
figure(2)
pcolor(log10(abs(real(plottestv)))), shading flat
colorbar