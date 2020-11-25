%% Processing script related to LS_2D_rt.m

clear all

outdir = 'output/test_nl_xy_01/';

load([outdir 'misc'])

intert = 1; % processing one every intert files
iteini = 1; % processing starts here
itefin = nt; % processing stops here

mass = zeros(1,nt/ppskip/intert+1);
angmom = mass;
energy = mass;
time = 0:dt*intert*ppskip:T;

for pp = iteini-1:ppskip*intert:itefin

    if size(num2str(pp),2) == np
        load([outdir 'psi_' num2str(pp)]);
    elseif size(num2str(pp),2) == np-1
        load([outdir 'psi_0' num2str(pp)]);
    elseif size(num2str(pp),2) == np-2
        load([outdir 'psi_00' num2str(pp)]);
    elseif size(num2str(pp),2) == np-3
        load([outdir 'psi_000' num2str(pp)]);
    elseif size(num2str(pp),2) == np-4
        load([outdir 'psi_0000' num2str(pp)]);
    elseif size(num2str(pp),2) == np-5
        load([outdir 'psi_00000' num2str(pp)]);
    end

    figure(21)
    subplot(1,4,[1 3])
    pcolor(X,Y,abs(wf).^2), shading flat
    title(['|\psi|^2, t = ' num2str(t)])
    axis equal tight
    colorbar
%    caxis([0 ])
    % hold off

    % [gradwft,gradwfr] = gradient(wf,dth,dr); % warning: not the gradient in polar coord.

    mass(pp/ppskip/intert+1) = sum(sum(abs(wf).^2))*dx*dy;
    %angmom(pp/ppskip/intert+1) = -i*sum(sum(conj(wf).*gradwft))*dr*dth;
    %energy(pp/ppskip/intert+1) = sum(sum(Rad(2:nr,:).*abs(gradwfr(2:nr,:)).^2 + ...
    %    (abs(gradwft(2:nr,:)).^2)./Rad(2:nr,:) + 0.5*abs(wf(2:nr,:)).^4))*dr*dth;

    %figure(21)
    subplot(1,4,4)
    plot(time,mass,'+'), shading flat
    axis tight
    title(['time evolution of the total mass'])

%     print('-f21','-djpeg','-r75',[outdir 'pics/pic_' num2str(pp)])

end

figure(22)
plot(0:dt*ppskip*intert:T,mass)
title('total mass')
% figure(23)
% plot(0:dt*ppskip*intert:T,angmom)
% title('total angular momentum')
% figure(24)
% plot(0:dt*ppskip*intert:T,energy)
% title('total energy')
