%% Processing script related to NLS_2D_rt.m

clear all

outdir = 'output/boundary/256_128_02/';
if ~exist([outdir 'pics'],'dir')
    mkdir([outdir 'pics'])
end
ftsize = 14;

load([outdir 'output/misc'])
%outdir = '../output/boundary/256_128_02/';
%ppskip = 1;

Theti = [Thet Thet(:,1)];
Radi = [Rad Rad(:,1)];

X = Radi.*cos(Theti);
Y = Radi.*sin(Theti);

intert = 4; % processing one every intert files
iteini = 0; % processing starts here
itefin = 25000; % processing stops here

%np = 3; % # of digits indexing the 'psi' files
mass = zeros(1,(itefin-iteini)/ppskip/intert+1);
angmom = mass;
energy = mass;
time = dt*(iteini:ppskip*intert:itefin);

for pp = iteini:ppskip*intert:itefin

    if size(num2str(pp),2) == np
        load([outdir 'output/psi_' num2str(pp)]);
    elseif size(num2str(pp),2) == np-1
        load([outdir 'output/psi_0' num2str(pp)]);
    elseif size(num2str(pp),2) == np-2
        load([outdir 'output/psi_00' num2str(pp)]);
    elseif size(num2str(pp),2) == np-3
        load([outdir 'output/psi_000' num2str(pp)]);
    elseif size(num2str(pp),2) == np-4
        load([outdir 'output/psi_0000' num2str(pp)]);
    elseif size(num2str(pp),2) == np-5
        load([outdir 'output/psi_00000' num2str(pp)]);
    elseif size(num2str(pp),2) == np-6
        load([outdir 'output/psi_000000' num2str(pp)]);
    elseif size(num2str(pp),2) == np-7
        load([outdir 'output/psi_0000000' num2str(pp)]);
    end

    wfi = [wf wf(:,1)];

    figure(21)
    %subplot(3,4,[1:3 5:7 9:11])
    subplot(4,3,[1:3 4:6 7:9])
    colormap(jet(256))
    surf(X,Y,abs(wfi).^2), shading interp
    xlabel('X','FontSize',ftsize)
    ylabel('Y','FontSize',ftsize)
    title(['2|\psi|^2, t = ' num2str(t)],'FontSize',ftsize)
    axis equal tight
    %caxis([0 1e-2])
    colorbar
    hold off

    [gradwft,gradwfr] = gradient(wf,dth,dr); % warning: not the gradient in polar coord.

    mass(pp/ppskip/intert+1) = sum(sum(Rad.*abs(wf).*abs(wf)))*dr*dth;
    angmom(pp/ppskip/intert+1) = -i*sum(sum(Rad.*conj(wf).*gradwft))*dr*dth;
    energy(pp/ppskip/intert+1) = dr*dth*sum(sum(...
        0.5*Rad(2:nr,:).*abs(gradwfr(2:nr,:)).^2 ...
        + 0.5*(abs(gradwft(2:nr,:)).^2)./Rad(2:nr,:) ...
        + 0.25*Rad(2:nr,:).*abs(wf(2:nr,:)).^4 ...
        + V*Rad(2:nr,:).*abs(wf(2:nr,:)).^2));

    %subplot(3,4,4)
    subplot(4,3,10)
    plot(time,mass,'+')
    ylabel('N','FontSize',ftsize,'Rotation',90)
    xlabel('time','FontSize',ftsize)
    axis([min(time) max(time) 0 500])
    %subplot(3,4,12)
    subplot(4,3,11)
    plot(time,real(angmom),'+')
    ylabel('L','FontSize',ftsize,'Rotation',90)
    xlabel('time','FontSize',ftsize)
    axis([min(time) max(time) 0 140])
    %subplot(3,4,8)
    subplot(4,3,12)
    plot(time,energy,'+')
    axis([min(time) max(time) 0 3.5])
    ylabel('E','FontSize',ftsize,'Rotation',90)
    xlabel('time','FontSize',ftsize)

    if pp == 0
        print('-f21','-depsc',[outdir 'pics/pic_' num2str(pp) '.eps'])
    end
    print('-f21','-djpeg',[outdir 'pics/pic_' num2str(pp) '.jpg'])

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
