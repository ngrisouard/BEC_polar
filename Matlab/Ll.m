function mat = Ll(n,r)

R = [0 0 r];

mat = sparse((n*n/2)*diag(R));
%mat = (n*n/2)*diag(R);


mat(1,4) = (-1)^n;
mat(2,3) = (-1)^n;

