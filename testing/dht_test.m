clear all

load('dht.mat')

ntn = 100;
R = 10;
nt = ntn;
dt = 5e-3; T = nt*dt;
ppskip = 40; % interval between two written outputs

nth = 256; dth = 2*pi/nth; % angular points
theta = 0:dth:2*pi-dth;

nr = 256; dr = R/(nr-1); % radial points
r = 0:dr:R;
interp_sch = 'spline';

jmodes = 128; % number of modes

[Thet,Rad] = meshgrid(theta,r);


    r1 = 3*R/4; % center of vortex (radius)
    theta1 = pi; % center of vortex (angle)
    circ1 = 2; % vortex circulation

    %r2 = R/4; % center of vortex (radius)
    %theta2 = pi; % center of vortex (angle)
    %circ2 = 1; % vortex circulation

    %r3 = 3*R/8; % center of vortex (radius)
    %theta3 = 0; % center of vortex (angle)
    %circ3 = 1; % vortex circulation
    
    
    wf = tanh((R-Rad)/sqrt(2)).*exp(1i*0*Thet)... % background
        .*(Rad.*exp(1i*circ1*Thet) - r1*exp(1i*circ1*theta1))... % 1st vortex
        ./sqrt(Rad.^2 + r1.^2 - 2*r1*Rad.*cos(circ1*(Thet-theta1)) + 1);
%wf = tanh((R-Rad)/sqrt(2))%.*exp(1i*7*Thet);
%wf = besselj(0,Rad.*2.4048./R).*exp(1i*7*Thet);

[indi,indj] = find(isinf(wf));      %set infinite value to 0
wf(indi,indj) = 0;

Theti = [Thet Thet(:,1)];
Radi = [Rad Rad(:,1)];

X = Radi.*cos(Theti);
Y = Radi.*sin(Theti);

wfi = [wf wf(:,1)];

figure(1)
colormap(jet(256))
%Z = angle(wfi);
Z = abs(wfi).^2;
pcolor(X,Y,Z), shading interp
axis equal tight
xlabel('x')
ylabel('y')
colorbar


%% Compute kernel

    comp_ker = 0;
    if comp_ker == 1
    %order = [1:nth/2 1:nth/2]
    %order_index = 1
        for ii = -nth/2:nth/2
            [H,kk,rr,I,KK,RR]=dht_mod([],R,jmodes,ii);
            save(['output/ker_' num2str(floor(ii))],'I','kk','rr','KK','RR')
            %order_index = order_index+1
        end

    end


%% fft, dht, idht, ifft

wf_old = wf;
[kt,fr] = obngfft(Thet,wf,2);
fr_old = fr;
%x = []
    for ii = -nth/2+1:nth/2
        
        load(['output/ker_' num2str(round(ii))])
        
        if abs(ii)>0
            I(1, 1) = 1;
            %kk(1) = 1/jmodes;
            %KK(1) = 1/jmodes;
            %RR(1) = 1/jmodes;
        end
        %kk(1) = 1/jmodes;
        % forward Hankel transform
        Fr = interp1(r,fr(:,ii+nth/2),rr,interp_sch,'extrap');
        adht = dht(Fr,RR,KK,I);
        Fr_old = Fr;
        
        %x= [x, Fr_old(1)-Fr(1)];
        % inverse Hankel transform

        %Fr = idht(adht.*exp(-1i*0.5*kk.^2*dt),I,KK,RR);
        Fr = idht(adht,I,KK,RR);
        %Fr(1) = Fr_old
        fr(:,ii+nth/2) = interp1(rr,Fr,r,interp_sch,'extrap'); 
        
        %Fr_new = idht(adht,I,KK,RR);
        
        %fr_new(:,ii+nth/2) = interp1(rr,Fr_new,r,interp_sch,'extrap');
        

    end

%[r2,wf] = obngifft(kt,fr,2);

%fr_new = [fr(:,2:end), 0];
%fr(1,:) = fr(1,128);
[r2,wf] = obngifft(kt,fr,2);

wfi = [wf wf(:,1)];
Z = abs(wfi).^2;
%Z = angle(wfi);
figure(2)
colormap(jet(256))
pcolor(X,Y,Z), shading interp
axis equal tight
xlabel('x')
ylabel('y')
colorbar

%%
%{
clear all
R = 10;
jmodes = 128;
[H,kk,rr,I,KK,RR]=dht_mod([],R,jmodes,0);
%}
