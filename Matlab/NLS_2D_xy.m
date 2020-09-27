%% Routine to solve the NLS equation in 2D and polar coordinates


clear all
tic

presdir = pwd;
outdir = 'output/test_nl_xy_01/';

%cd(outdir)

split_NL = 1; % =0.5 for NLS, =1 for LS, =1/3 for NLS with dissipation
V = -0.5;

L = 10; % length of one size of the box NORMALIZED BY THE HEALING LENGTH

tic





%% Time and spatial coordinates

nti = 0; % starting iteration (loads the last file of the previous run)
nt = 1; dt = 1e-2; T = nt*dt;
np = length(num2str(nt));
ppskip = 1; % interval between two written outputs
if nti > 0
    comp_ker = 0;
end

ny = 256; dy = L/ny;
y = -L/2+dy:dy:L/2;
ky = k_of_x(y);

nx = 256; dx = L/nx;
x = -L/2+dx:dx:L/2;
kx = k_of_x(x);

%jmodes = nr;

[Y,X] = meshgrid(y,x);
[KY,KX] = meshgrid(ky,kx);




%% Initial condition

%wf = ones(nx,ny);
wf = (1-exp(-X.^2 -Y.^2)).*(X + i*Y)./sqrt(X.^2 + Y.^2);
wf(ny/2,nx/2) = 0;
%wf = cos(X/L)+i*sin(Y/L);








%% Temporal loop

pp = 1+nti; % indicator of the progression in the loop

tic


if nti > 0

    if size(num2str(nti),2) == np
        load([outdir 'psi_' num2str(nti)]);
    elseif size(num2str(pp),2) == np-1
        load([outdir 'psi_0' num2str(nti)]);
    elseif size(num2str(pp),2) == np-2
        load([outdir 'psi_00' num2str(nti)]);
    elseif size(num2str(pp),2) == np-3
        load([outdir 'psi_000' num2str(nti)]);
    elseif size(num2str(pp),2) == np-4
        load([outdir 'psi_0000' num2str(nti)]);
    else
        error('Loaded file will be wrong.')
    end
    
else

    t = 0;
    if nt >= 1 && nt < 10
        save([outdir 'psi_0'],'wf','t')
    elseif nt >= 10 && nt < 100
        save([outdir 'psi_00'],'wf','t')
    elseif nt >= 100 && nt < 1000
        save([outdir 'psi_000'],'wf','t')
    elseif nt >= 1000 && nt < 10000
        save([outdir 'psi_0000'],'wf','t')
    elseif nt >= 10000 && nt < 100000
        save([outdir 'psi_00000'],'wf','t')
    else
        error('Output names wont be nice-looking, complete this section.')
    end

end


save([outdir 'misc'],'x','y','T','nt','dt','ny','dy','nx','dx','L',...
    'ppskip','X','Y','np')

optime = exp(-i*0.5*split_NL*dt*(KX.^2+KY.^2)*0);

for t = (nti+1)*dt:dt:T
    
    % 1st step: the laplacian part
    figure(21)
    subplot(1,2,1)
    pcolor(abs(wf).^2),shading flat
    colorbar
    axis equal tight
    WF = fft2(wf);
    subplot(1,2,2)
    pcolor(abs(WF).^2),shading flat
    colorbar
    axis equal tight
    wf = ifft2(WF);

    % 2nd step: NL part
    if split_NL == 0.5 || split_NL == 0 || abs(split_NL-1/3) < 1e-13

        wf = wf.*exp(-i*(V + 0.5*abs(wf).^2)*dt*split_NL);

    end


     % 3rd step: dissipation

    if split_NL == 0 || abs(split_NL-1/3) < 1e-13

        error('Dissipation not implemented yet.');

    end

    if abs(pp/ppskip - floor(pp/ppskip)) <= 1e-14
        % in this case, an output file is produced
        % name of the outputfile
        if size(num2str(pp),2) == np
            nameout = [outdir 'psi_' num2str(pp)];
        elseif size(num2str(pp),2) == np-1
            nameout = [outdir 'psi_0' num2str(pp)];
        elseif size(num2str(pp),2) == np-2
            nameout = [outdir 'psi_00' num2str(pp)];
        elseif size(num2str(pp),2) == np-3
            nameout = [outdir 'psi_000' num2str(pp)];
        elseif size(num2str(pp),2) == np-4
            nameout = [outdir 'psi_0000' num2str(pp)];
        else
            error('Output names will be wrong.')
        end
        save(nameout,'wf','t')
    end
    pp = pp+1;

end

toc

beep

%cd(al)