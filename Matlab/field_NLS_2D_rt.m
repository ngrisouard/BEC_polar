%% Processing script related to NLS_2D_rt.m

clear all

outdir = 'output/valid/capillary/output/';
ftsize = 28;

load([outdir 'misc'])
%outdir = '../output/boundary/256_128_02/';
%ppskip = 1;

Theti = [Thet Thet(:,1)];
Radi = [Rad Rad(:,1)];

X = Radi.*cos(Theti);
Y = Radi.*sin(Theti);

intert = 1; % processing one every intert files
iteini = 60; % processing starts here
itefin = 60; % processing stops here

% np = 3; % # of digits indexing the 'psi' files
mass = zeros(1,(itefin-iteini)/ppskip/intert+1);
angmom = mass;
energy = mass;
warning('Check if the time is right or not.')
time = dt*(iteini:ppskip*intert:itefin);

for pp = iteini:ppskip*intert:itefin

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
    elseif size(num2str(pp),2) == np-6
        load([outdir 'psi_000000' num2str(pp)]);
    elseif size(num2str(pp),2) == np-7
        load([outdir 'psi_0000000' num2str(pp)]);
    end

    wfi = [wf wf(:,1)];

    figure(21)
    % subplot(3,4,[1:3 5:7 9:11])% ; 5 7 ; 9 11])
    colormap(hot(256))
    set(gca,'FontSize',ftsize);
    %set(gca,'DataAspectRatio',[1 1 1])
    %set(h1,'ZTick',[0 1],'ZTicklabel',{'0' '1'});
    surf(X,Y,abs(wfi).^2), shading interp
    daspect([1 1 0.75])
    axis([-10 10 -10 10 0 1.6])
    xlabel('X','FontSize',ftsize)
    ylabel('Y','FontSize',ftsize)
    title(['|\psi|^2, t = ' num2str(t)],'FontSize',ftsize)
    %axis equal tight
    caxis([0.5 1.4])

    %colorbar
    hold off
%
%     [gradwft,gradwfr] = gradient(wf,dth,dr); % warning: not the gradient in polar coord.
%
%     mass(pp/ppskip/intert+1) = sum(sum(Rad.*abs(wf).*abs(wf)))*dr*dth;
%     angmom(pp/ppskip/intert+1) = -i*sum(sum(Rad.*conj(wf).*gradwft))*dr*dth;
%     energy(pp/ppskip/intert+1) = dr*dth*sum(sum(...
%         0.5*Rad(2:nr,:).*abs(gradwfr(2:nr,:)).^2 ...
%         + 0.5*(abs(gradwft(2:nr,:)).^2)./Rad(2:nr,:) ...
%         + 0.25*Rad(2:nr,:).*abs(wf(2:nr,:)).^4 ...
%         + V*Rad(2:nr,:).*abs(wf(2:nr,:)).^2));
%     % term1(pp/ppskip/intert+1) = sum(sum(0.5*Rad(2:nr,:).*abs(gradwfr(2:nr,:)).^2));
%     % term2(pp/ppskip/intert+1) = sum(sum(0.5*(abs(gradwft(2:nr,:)).^2)./Rad(2:nr,:)));
%     % term3(pp/ppskip/intert+1) = sum(sum(0.25*Rad(2:nr,:).*abs(wf(2:nr,:)).^4));
%     % term4(pp/ppskip/intert+1) = sum(sum(V*Rad(2:nr,:).*abs(wf(2:nr,:)).^2));
%
%     subplot(3,4,4)
%     plot(time,mass,'+')%,'markersize',1), shading flat
%     title('total mass')
%     xlabel('time')
%     subplot(3,4,8)
%     plot(time,real(angmom),'+')%,'markersize',1)
%     title('total angular momentum')
%     xlabel('time')
%     subplot(3,4,12)
%     plot(time,energy,'+')%,'markersize',1)
%     %plot(time,term1,'b')%,...
%         %time,term3,'r')
%         %'.','markersize',1)
%     title('total energy')
%     xlabel('time')
%
% % term4 = mass;
%    M(pp/ppskip/intert+1) = getframe;

    %print('-f21','-djpeg',[outdir 'pics/pic_' num2str(pp)])
    print('-f21','-depsc','../figures/capillaryBEC.eps')

end

% movie(M)

% figure(22)
% plot(time,mass)
% title('total mass')
% figure(23)
% plot(time,angmom)
% title('total angular momentum')
% figure(24)
% plot(time,energy)
% title('total energy')
