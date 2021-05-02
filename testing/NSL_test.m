nti = 0; % starting iteration (loads the last file of the previous run)
ntn = 100; % final iteration of the current run

if nti == 0 % then define all parameters

    comp_ker = 1; % =1 if the values of the kernel have to be computed

    split_NL = 0.5; % =0.5 for NLS, =1 for LS, =1/3 for NLS with dissipation
    if split_NL == 1
        V = 0;
    else
        V = -0.5;
    end

    R = 10; % radius of the cup NORMALIZED BY THE HEALING LENGTH



    %% Time and spatial coordinates

    nt = ntn;
    dt = 5e-3; T = nt*dt;
    ppskip = 40; % interval between two written outputs
    if nti > 0
        comp_ker = 0;
    end

    nth = 256; dth = 2*pi/nth; % angular points
    theta = 0:dth:2*pi-dth;

    nr = 256; dr = R/(nr-1); % radial points
    r = 0:dr:R;
    interp_sch = 'cubic';

    jmodes = 128; % number of modes

    [Thet,Rad] = meshgrid(theta,r);



    %% Initial condition

    r1 = 0; % center of vortex (radius)
    theta1 = 0; % center of vortex (angle)
    circ1 = 1; % vortex circulation

    r2 = 3*R/8; % center of vortex (radius)
    theta2 = pi; % center of vortex (angle)
    circ2 = -1; % vortex circulation

    %r3 = 3*R/8; % center of vortex (radius)
    %theta3 = 0; % center of vortex (angle)
    %circ3 = 1; % vortex circulation


    wf = tanh((R-Rad)/sqrt(2)).*exp(1i*7*Thet)... % background
        .*(Rad.*exp(1i*circ1*Thet) - r1*exp(1i*circ1*theta1))... % 1st vortex
        ./sqrt(Rad.^2 + r1.^2 - 2*r1*Rad.*cos(circ1*(Thet-theta1)) + 1);
    if exist('r2','var')
        wf = wf.*(Rad.*exp(1i*circ2*Thet) - r2*exp(1i*circ2*theta2))... % 2nd vortex
            ./sqrt(Rad.^2 + r2.^2 - 2*r2*Rad.*cos(circ2*(Thet-theta2)) + 1);
    end
    if exist('r3','var')
        wf = wf.*(Rad.*exp(1i*Thet) - r3*exp(1i*theta3))... % 3rd vortex
            ./sqrt(Rad.^2 + r3.^2 - 2*r3*Rad.*cos(circ3*(Thet-theta3)) + 1);
    end

    [indi,indj] = find(isinf(wf));      %set infinite value to 0
    wf(indi,indj) = 0;

    Theti = [Thet Thet(:,1)];
    Radi = [Rad Rad(:,1)];

    X = Radi.*cos(Theti);
    Y = Radi.*sin(Theti);

    wfi = [wf wf(:,1)];

    figure(1)
    colormap(jet(256))
    surf(X,Y,abs(wfi).^2), shading interp
    axis equal tight
    xlabel('x')
    ylabel('y')
    colorbar
end


pp = 1+nti; % indicator of the progression in the loop


deriv_thet = 1i*(ones(nr,1)*(-nth/2+1:nth/2));

[kt,fr] = obngfft(Thet,wf,2); % fft to isolate every angular mode
    if abs((pp-1)/ppskip - floor((pp-1)/ppskip)) <= 1e-14
        % derivative of wf % theta FOR THE PREVIOUS ITERATION
        [r2,dwfdt] = obngifft(kt,deriv_thet.*fr,2);
    end
    
    ii = -10
    load(['output/ker_' num2str(round(ii))])

        % forward Hankel transform
    Fr = interp1(r,fr(:,ii+nth/2),rr,interp_sch,'extrap');
        adht = dht(Fr,RR,KK,I);

        % inverse Hankel transform
        Fr = idht(adht.*exp(-1i*0.5*kk.^2*dt),I,KK,RR);
        fr(:,ii+nth/2) = interp1(rr,Fr,r,interp_sch,'extrap');