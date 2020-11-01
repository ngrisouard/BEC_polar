clear all

R = 10
jmodes = 128
ii = 215
N=256
n = ii
h=[]
%{
if isempty(N)
   N=256;
end
if isempty(n)
   n=0;
end
if numel(n) > 1
   K=N;
   I=n;
else
   if ~isempty(h) & isa(h,'numeric')
      error('Need a function h(r) without kernel.');
   end
   load('dht.mat');                 % Bessel Jn rooths
   C=c(1+n,1+N);
   c=c(1+n,1:N);
   r=R/C*c(:);
   k=c(:)/R;
   I=abs(besselj(1+n,c));
   if n > 0
      I(1)=1/N;                     % avoid zero - thanks to Nicolas Grisouard
   end
%    I(~I)=1/N;                     % or this, but there should be only one
   K=2*pi*R/C*I(:);
   R=I(:)/R;
   I=sqrt(2/C)./I;
   I=I(:)*I.*besselj(n,c(:)/C*c);
end
if isempty(h)
   H=h;
else
   if ~isa(h,'numeric')
      h=feval(h,r);
   end
   H=I*(h./R).*K; 
end
%}

[H,k,r,I,K,R, h]=dht([],R,jmodes,ii)