clear all

%% Program to solve the NLS in 2D
% psi_t = (i/2)*lapl(psi) - i*K*(f(r)+|psi|^2)*psi
tic

outdir = '../output/test_pw_x_01/';
order = 1;

xi = 1e-2; % healing length
n0 = 1; % background density
nl = 0; % 0 to shut down the nonlinear term

%% physical domain
Lx = 2.; Ly = Lx; % length of one side of the domain
nx = 256; dx = Lx/nx;
ny = 8; dy = Ly/ny;
x = -Lx/2+dx:dx:Lx/2;
y = -Ly/2+dy:dy:Ly/2;
[X,Y] = meshgrid(x,y);

ix = 1:nx; iy = 1:ny;

%r = sqrt(x.^2 + y.^2);
f = 1e8*ones(ny,nx)*0;
% for k=1:ny
%     for j=1:nx
%         if sqrt(X(k,j)^2 + Y(k,j)^2) <= 1
%             f(k,j) = 0; % potential
%         end
%     end
% end
%f = ones(nx,ny);
%f = 1./f1 - 1; % potential
%figure(3)
%pcolor(X,Y,f), shading flat
%colorbar
%axis equal tight

k_x = pi/Lx; k_y = pi/Ly;
%psi_0 = cos(k_x*X);
psi_0 = exp(-((1+i)*(X/0.4).^2));
psi_0 = psi_0/(Lx*Ly*sum(sum(psi_0.*conj(psi_0))));
omega = k_x^2/4;

T = k_x^2; % duration
dt = k_x^2/4000; % time step

figure(1)
pcolor(X,Y,psi_0.*conj(psi_0)), shading flat
axis equal tight


%% iterations
psi2_time = zeros(1,T/dt+1);
time = zeros(1,T/dt+1);
tt = 0;
N = -i*psi_0.*conj(psi_0)/(2*xi^2*n0)*nl;

if order == 2
    
    [kx,ky,psi] = obfft2(x,y,psi_0);
    [Kx,Ky] = meshgrid(kx,ky);
    D = i*(Kx.^2 + Ky.^2)/2 - i*f/(2*xi^2*n0);
    [x,y,psi] = obifft2(kx,ky,exp(dt*D/2).*psi);
    [kx,ky,psi] = obfft2(x,y,exp(N*dt).*psi);
    [x,y,psi] = obifft2(kx,ky,exp(dt*D/2).*psi);
    psi2_time(1) = psi_0(ny/2,nx/2)*conj(psi_0(ny/2,nx/2));

    for tt = 1:T/dt;
        N = -i*psi.*conj(psi)/(2*xi^2*n0)*nl;
        [kx,ky,psi] = obfft2(x,y,psi);
        [Kx,Ky] = meshgrid(kx,ky);
        [x,y,psi] = obifft2(kx,ky,exp(dt*D/2).*psi);
        [kx,ky,psi] = obfft2(x,y,exp(N*dt).*psi);
        [x,y,psi] = obifft2(kx,ky,exp(dt*D/2).*psi);
        psi2_time(tt+1) = psi(ny/2,nx/2)*conj(psi(ny/2,nx/2));
        time(tt+1) = tt*dt;
    end

elseif order == 1
    
    [kx,ky,psi] = obfft2(x,y,exp(dt*N).*psi_0);
    [Kx,Ky] = meshgrid(kx,ky);
    D = i*(Kx.^2 + Ky.^2)/2 - i*f/(2*xi^2*n0);
    [x,y,psi] = obifft2(kx,ky,exp(dt*D).*psi);
    psi2_time(1) = psi_0(ny/2,nx/2)*conj(psi_0(ny/2,nx/2));

    for tt = 1:T/dt;
        N = -i*psi.*conj(psi)/(2*xi^2*n0)*nl;
        [kx,ky,psi] = obfft2(x,y,psi.*exp(dt*N));
        [x,y,psi] = obifft2(kx,ky,exp(dt*D).*psi);
        psi2_time(tt+1) = psi(ny/2,nx/2)*conj(psi(ny/2,nx/2));
        time(tt+1) = tt*dt;
    end

end

figure(2)
plot(time,psi2_time)
%pcolor(X,Y,psi.*conj(psi)), shading flat
%axis equal tight

toc