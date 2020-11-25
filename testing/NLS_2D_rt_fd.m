%% Routine to solve the NLS equation in 2D and polar coordinates


clear all
tic

outdir = '../output/256_256_128_01/';

nti = 0; % starting iteration (loads the last file of the previous run)
ntn = 100; % final iteration of the current run

if nti == 0 % then define all parameters

    split_NL = 0.5; % =0.5 for NLS, =1 for LS, =1/3 for NLS with dissipation
    V = -0.5;

    R = 10; % radius of the cup NORMALIZED BY THE HEALING LENGTH



    %% Time and spatial coordinates

    nt = ntn; dt = 1e-2; T = nt*dt;
    ppskip = 10; % interval between two written outputs
    if nti > 0
        comp_ker = 0;
    end

    nth = 16; dth = 2*pi/nth;
    theta = 0:dth:2*pi-dth;

    nr = 16; dr = 2*R/(2*nr+1);
    r = dr*(1:(nr+1) - 0.5);

    [Thet,Rad] = meshgrid(theta,r);



    %% Initial condition

    % % scaling of the radial modes
    % zero_bess_t = zeros(10,nth);
    %
    % zero_bess_t(:,1) = bessel_zeros(1,0,10,1e-30);
    % zero_bess_t(1,nth/2+1:nth) = 0;
    % for n = 1:nth/2
    %     zero_bess_t(2:10,n+nth/2) = bessel_zeros(1,n,9,1e-30);
    % end
    % for n = 1:nth/2-1
    %     zero_bess_t(:,n) = zero_bess_t(:,nth-n);
    % end
    % K = zero_bess_t/R;
    %
    % m = 1;
    % l = 5;
    % wf = besselj(m,K(l,m+nth/2)*Rad).*exp(i*m*Thet)...
    %         +besselj(m-5,K(l+5,m-5+nth/2)*Rad).*exp(i*(m-5)*Thet);
    %         +besselj(m+5,K(l+5,m+5+nth/2)*Rad).*exp(i*(m+5)*Thet)...
    % wf = Rad.*(R-Rad).*exp(i*Thet);
    % wf = wf/sqrt(sum(sum(Rad.*abs(wf).^2))*dr*dth);
    % wf = Rad - R/2 +i*(Thet - pi);
    % wf = ones(nr,nth);
    % wf = tanh((R-Rad)/sqrt(2));
    % wf(1,:) = 0;

    r0 = 3*R/4;
    theta0 = pi;
    wf = tanh((R-Rad)/sqrt(2))...
        .*tanh(sqrt((Rad.^2 + r0.^2 - 2*r0*Rad.*cos(Thet-theta0))/2))...
        .*(Rad.*exp(i*Thet) - r0*exp(i*theta0))...
        ./sqrt((Rad.^2 + r0.^2 - 2*r0*Rad.*cos(Thet-theta0))/2);

    ind = find(isinf(wf));
    wf(ind) = 0;


    M1 = sparse(zeros(nr+2));
    %M1 = zeros(nr+2);
    M2 = M1;
    for ii = 3:nr

        M1(ii,ii) = 30;

        M1(ii,ii-1) = -16;
        M2(ii,ii-1) = 8;

        M1(ii,ii+1) = -16;
        M2(ii,ii+1) = -8;

        M1(ii,ii-2) = 1;
        M2(ii,ii-2) = -1;

        M1(ii,ii+2) = 1;
        M2(ii,ii+2) = 1;
    end

    M1(nr+1,nr-2) = 1;
    M2(nr+1,nr-2) = 1;

    M1(nr+1,nr-1) = -4;
    M2(nr+1,nr-1) = -6;

    M1(nr+1,nr) = -6;
    M2(nr+1,nr) = 18;

    M1(nr+1,nr+1) = 20;
    M2(nr+1,nr+1) = -10;

    M1(nr+1,nr+2) = -11;
    M2(nr+1,nr+2) = -3;

    RM = sparse([ones(2,nr+2) ; (1./r(:))*ones(1,nr+2)]);
    %RM = [ones(2,nr+2) ; (1./r(:))*ones(1,nr+2)];

    M1 = M1/(24*dr*dr);
    M2 = M2.*RM/(24*dr);

    M = M1+M2;

    iI = zeros(nr+2);
    iI(3:nr+2,3:nr+2) = i*eye(nr);
    %iI = sparse(iI);

    t = 0;
    save([outdir 'psi_00000000'],'wf','t')
    np = 8;

    % In this file, all the parameters are saved, not just those which will be
    % used later. It can be useful if there is a doubt on the parameters used.
    save([outdir 'misc'],'r','theta','T','nt','dt','nth','dth','nr','dr',...
        'R','ppskip','Rad','Thet','np','V','outdir','wf','t','split_NL'...
        ,'M','iI')


else % then reload all the parameters and last iteration

    load([outdir 'misc.mat'])
    nt = ntn;
    T = dt*nt;
    save([outdir 'misc'],'r','theta','T','nt','dt','nth','dth','nr','dr',...
        'R','ppskip','Rad','Thet','np','V','outdir','split_NL','wf','t',...
        'M','iI')

    if size(num2str(nti),2) == np
        load([outdir 'psi_' num2str(nti)]);
    elseif size(num2str(nti),2) == np-1
        load([outdir 'psi_0' num2str(nti)]);
    elseif size(num2str(nti),2) == np-2
        load([outdir 'psi_00' num2str(nti)]);
    elseif size(num2str(nti),2) == np-3
        load([outdir 'psi_000' num2str(nti)]);
    elseif size(num2str(nti),2) == np-4
        load([outdir 'psi_0000' num2str(nti)]);
    elseif size(num2str(nti),2) == np-5
        load([outdir 'psi_00000' num2str(nti)]);
    elseif size(num2str(nti),2) == np-6
        load([outdir 'psi_000000' num2str(nti)]);
    elseif size(num2str(nti),2) == np-7
        load([outdir 'psi_0000000' num2str(nti)]);
    else
        error('Loaded file will be wrong.')
    end

end



%% Temporal loop

pp = 1+nti; % indicator of the progression in the loop

tic


if nti > 0

end

for t = (nti+1)*dt:dt:T

    %t
    %mass = sum(sum(abs(wf).^2.*Rad))*dr*dth

    % 1st step: the linear step, computation of the laplacian part.

    [kt,fr] = obngfft(Thet,wf,2); % fft to isolate every angular mode
    fr_int = zeros(nth,nr+2);
    fr_iint = (iI + 0.5*M)*fr_int';

    %massr = sum(sum(abs(fr).^2.*Rad))*dr/2/pi

    for ii = -nth/2+1:nth/2
        
        fr_int(ii+nth/2,:) = [fr(2,ii+nth/2)*(-1)^ii ;...
            fr(1,ii+nth/2)*(-1)^ii ;...
            fr(:,ii+nth/2)];
        fr_int(ii+nth/2,:) = (iI - 0.5*(M+Ll(ii,r))) \...
            (fr_iint(:,ii+nth/2) + 0.5*Ll(ii,r)*fr_int(ii+nth/2,:)');
%         fr_int(ii+nth/2,:) = linsolve(iI - 0.5*(M+Ll(ii,r)),...
%             fr_iint(:,ii+nth/2) + 0.5*Ll(ii,r)*fr_int(ii+nth/2,:)');

    end
    
    fr = fr_int(:,3:nr+2)';

    %massrrrr = sum(sum(abs(fr).^2.*Rad))*dr/2/pi

    [r2,wf] = obngifft(kt,fr,2);

    % 2nd step: NL part

    if split_NL == 0.5 || split_NL == 0 || abs(split_NL-1/3) < 1e-13

        %mass_L = sum(sum(abs(wf).^2.*Rad))*dr*dth
        wf = wf.*exp(-i*(V + 0.5*abs(wf).^2)*dt*split_NL);
        %mass_NL = sum(sum(abs(wf).^2.*Rad))*dr*dth

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
        elseif size(num2str(pp),2) == np-5
            nameout = [outdir 'psi_00000' num2str(pp)];
        elseif size(num2str(pp),2) == np-6
            nameout = [outdir 'psi_000000' num2str(pp)];
        elseif size(num2str(pp),2) == np-7
            nameout = [outdir 'psi_0000000' num2str(pp)];
        else
            error('Output names will be wrong.')
        end
        save(nameout,'wf','t')
    end
    pp = pp+1;

end

toc

beep