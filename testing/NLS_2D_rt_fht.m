%% Routine to solve the NLS equation in 2D and polar coordinates

clear all
tic

outdir = 'output/test_15/';

comp_ker = 0; % =1 if the values of the kernel or the time operator have to be computed

split_NL = 1; % =0.5 for NLS, =1 for LS, =1/3 for NLS with dissipation
V = -1/2; %internal potential

R = 10; % radius of the cup NORMALIZED BY THE HEALING LENGTH

tic





%% Time and spatial coordinates

nti = 0; % starting iteration (loads the last file of the previous run)
nt = 100; dt = 1e-1; T = nt*dt;
ntot = nti+nt;
np = length(num2str(nt));
ppskip = 1; % interval between two written outputs
if nti > 0
    comp_ker = 0;
end

nth = 256; dth = 2*pi/nth;
theta = 0:dth:2*pi-dth;

nr = 256; dr = R/(nr-1);
r = 0:dr:R;
interp_sch = 'spline';

jmodes = nr;

[Thet,Rad] = meshgrid(theta,r);






%% Initial condition


m = 1;
l = 5;
% wf = besselj(m,K(l,m+nth/2)*Rad).*exp(i*m*Thet)...
%         +besselj(m-5,K(l+5,m-5+nth/2)*Rad).*exp(i*(m-5)*Thet);
%         +besselj(m+5,K(l+5,m+5+nth/2)*Rad).*exp(i*(m+5)*Thet)...
% elseif strcmp(outdir,'output/test_rt_02/') || strcmp(outdir,'output/test_nl_rt_02/')
%     wf = Rad.*(R-Rad).*exp(i*m*Thet);
% elseif strcmp(outdir,'output/test_rt_03/')
%     wf = besselj(m,K(l,m+nth/2)*Rad).*exp(i*m*Thet)...
%         +besselj(m-25,K(l+5,m-25+nth/2)*Rad).*exp(i*(m-25)*Thet);
% end
wf = Rad.*(R-Rad).*exp(i*m*Thet);
% wf = wf/sqrt(sum(sum(Rad.*abs(wf).^2))*dr*dth);
% wf = Rad - R/2 +i*(Thet - pi);








%% Computation of the integration kernel (optionnal)
% When the quasi fast Hankel transform is performed, the computation of the
% so-called "integration kernel" is a long task but does not depend on the
% input function. It can (has to) be done outside the temporal loop.

if comp_ker == 1
    
    % scaling of the radial modes
    zero_bess = zeros(jmodes,nth);

    zero_bess(:,1) = bessel_zeros(1,0,jmodes,1e-30);
    zero_bess(1,nth/2+1:nth) = 0;
    for n = 1:nth/2
        zero_bess(2:jmodes,n+nth/2) = bessel_zeros(1,n,jmodes-1,1e-30);
    end
    for n = 1:nth/2-1
        zero_bess(:,n) = zero_bess(:,nth-n);
    end
    K = zero_bess/R;

    % computation of the kernel (see Marcel Leutenegger or Siegman, Opt. Lett., 1977)
    mm = 10; % minimum pts per cycle
    MM = 14; % maximum pts per cycle
    a = max(K,[],1)*R/2/pi;
    N = pow2(ceil(1 + log2(a.*mm.*log(a*MM))));
    a = 1./a/mm;
    ro = R*exp(-a.*N/2);
    ko = max(K,[],1).*exp(-a.*N/2);
    for ii = nth/2:nth % only half of the modes +1 as e^(i*n*t) and e^(-i*n*t) are eq.
        N1 = N(ii);
        I = exp(a(ii)*(0:N1-1));
        kk = ko(ii)*I(1:N1/2);                         % samplings
        N2 = numel(kk);
        rr = ro(ii)*I(1:N1/2);
        if ii == nth || ii == nth/2
            I = ifft(a(ii)*ko(ii)*ro(ii)*I.*besselj(ii-nth/2,ko(ii)*ro(ii)*I));
            optime = exp(-i*0.5*kk.^2*dt*split_NL);
            save([outdir 'ker_' num2str(floor(ii-nth/2))],...
                'I','kk','rr','N1','N2','optime')
        else
            I = [ifft(a(ii)*ko(ii)*ro(ii)*I.*besselj(ii-nth/2,ko(ii)*ro(ii)*I));...
                 ifft(a(ii)*ko(ii)*ro(ii)*I.*besselj(-ii+nth/2,ko(ii)*ro(ii)*I))]; % kernel
            RR = ones(2,1)*rr;
            KK = ones(2,1)*kk;
            optime = exp(-i*0.5*KK.^2*dt*split_NL);
            save([outdir 'ker_' num2str(floor(ii-nth/2))],...
                'I','kk','rr','RR','KK','N1','N2','optime')
        end
    end
end






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
    if nt >= 10 && nt < 100
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

for t = (nti+1)*dt:dt:T
    
    massrr = 0;
    massk = 0;
    masskk = 0;
    massrrr = 0;

    % 1st step: the linear step, computation of the laplacian part.
    
    %t
    %mass = sum(sum(abs(wf).^2.*Rad))*dr*dth
    
    [kt,fr] = obngfft(Thet,wf,2); % fft to isolate every angular mode

    % Hankel transform is easier when performed along r for each mode =>
    % loop. It is possible to deal with two modes at a time tough, except
    % for the 0 and last modes due to the asymmetry of these two cases.

    load([outdir 'ker_0']);

    %massr = sum(sum(abs(fr).^2.*Rad))*dr/2/pi;
    Fr = interp1(r,fr(:,nth/2),rr,interp_sch);
    %massrr = massrr + sum(abs(Fr(2:N2)).^2.*rr(2:N2).*diff(rr))/2/pi;
    afht = fft(fft(Fr.*rr,N1).*I);           % transform
    if isreal(Fr)
        afht=real(afht);
    end
    afht = 2*pi*afht(1:N1/2)./kk;
    %massk = massk + sum(abs(afht(2:end)).^2.*kk(2:end).*diff(kk))/2/pi;
    %masskk = masskk + sum(abs(afht(2:end).*optime(2:end)).^2.*kk(2:end).*diff(kk))/2/pi;

    Fr = fft(fft(afht.*kk.*optime,2*N2).*I);        % inverse transform
    if isreal(afht.*optime)
        Fr=real(Fr);
    end
    Fr = Fr(1:N2)./(2*pi*rr);
    %massrrr = massrrr + sum(abs(Fr(2:N2)).^2.*rr(2:N2).*diff(rr))/2/pi;
    fr(:,nth/2) = interp1(rr,Fr,r,interp_sch);
    
    clear Fr afht rr kk
    
    for ii = 2:nth/2
        
        load([outdir 'ker_' num2str(floor(ii-1))])
        
        Fr = [interp1(r,fr(:,ii+nth/2-1),rr,interp_sch);...
            interp1(r,fr(:,-ii+nth/2+1),rr,interp_sch)];
        %massrr = massrr + sum(sum(abs(Fr(:,2:N2)).^2.*RR(:,2:N2).*diff(RR,[],2)))/2/pi;
        afht = fft(fft(Fr.*RR,N1,2).*I,[],2);     % transform
        if isreal(Fr)
            afht=real(afht);
        end
        afht = 2*pi*afht(:,1:N1/2)./KK;
        %massk = massk + sum(sum(abs(afht(:,2:end)).^2.*KK(:,2:end).*diff(KK,[],2)))/2/pi;
        %masskk = masskk + sum(sum(abs(afht(:,2:end).*optime(:,2:end)).^2.*KK(:,2:end).*diff(KK,[],2)))/2/pi;

        Fr = fft(fft(afht.*KK.*optime,2*N2,2).*I,[],2);
        %Fr = fft(fft(afht.*KK,2*N2,2).*I,[],2);
        if isreal(afht.*optime)
            Fr=real(Fr);
        end
        Fr = Fr(:,1:N2)./(2*pi*RR);
        %massrrr = massrrr + sum(sum(abs(Fr(:,2:N2)).^2.*RR(:,2:N2).*diff(RR,[],2)))/2/pi;
        fr(:,ii+nth/2-1) = interp1(rr,Fr(1,:),r,interp_sch);
        fr(:,-ii+nth/2+1) = interp1(rr,Fr(2,:),r,interp_sch);

        clear Fr afht rr kk RR KK

    end

    load([outdir 'ker_' num2str(floor(nth/2))]);

    Fr = interp1(r,fr(:,nth),rr,interp_sch);
    %massrr = massrr + sum(abs(Fr(2:N2)).^2.*rr(2:N2).*diff(rr))/2/pi
    afht = fft(fft(Fr.*rr,N1).*I);           % transform
    if isreal(Fr)
      afht=real(afht);
    end
    afht = 2*pi*afht(1:N1/2)./kk;
    %massk = massk + sum(abs(afht(2:end)).^2.*kk(2:end).*diff(kk))/2/pi
    %masskk = masskk + sum(abs(afht(2:end).*optime(2:end)).^2.*kk(2:end).*diff(kk))/2/pi

    Fr = fft(fft(afht.*kk.*optime,2*N2).*I);
    %Fr = fft(fft(afht.*kk,2*N2).*I);
    if isreal(afht.*optime)
        Fr=real(Fr);
    end
    Fr = Fr(1:N2)./(2*pi*rr);
    %massrrr = massrrr + sum(abs(Fr(2:N2)).^2.*rr(2:N2).*diff(rr))/2/pi
    fr(:,nth) = interp1(rr,Fr,r,interp_sch);
    
    %massrrrr = sum(sum(abs(fr).^2.*Rad))*dr/2/pi

    clear Fr afht rr kk

    [r2,wf] = obngifft(kt,fr,2);
    
    
    if split_NL == 0.5 || abs(split_NL-1/3) < 1e-13 % 2nd step: NL part
        
        mass_L = sum(sum(abs(wf).^2.*Rad))*dr*dth;
        wf = wf.*exp(-i*(V + 0.5*abs(wf).^2)*dt*split_NL);
        mass_NL = sum(sum(abs(wf).^2.*Rad))*dr*dth;
        
    end
    
    if abs(split_NL-1/3) < 1e-13 % 3rd step: dissipation
        
        error('Dissipation not implemented yet.')
        
    end
    
    if abs(pp/ppskip - floor(pp/ppskip)) < 1e-13
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

save([outdir 'misc'],'r','theta','T','nt','dt','nth','dth','nr','dr','R',...
    'ppskip','Rad','Thet','ntot')

toc

beep