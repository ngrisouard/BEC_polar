
clear all

if ~exist('output','dir')
    mkdir output
end

if ~exist('dht.m','file')
    !cp ~/../../Volumes/LaCie/Recherche/GFD/Project/code/dht.m .
    !cp ~/../../Volumes/LaCie/Recherche/GFD/Project/code/idht.m .
    !cp ~/../../Volumes/LaCie/Recherche/GFD/Project/code/obng*.m .
    !cp ~/../../Volumes/LaCie/Recherche/GFD/Project/code/x_of_k_ng.m .
    !cp ~/../../Volumes/LaCie/Recherche/GFD/Project/code/k_of_x_ng.m .
end

tic

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
    interp_sch = 'spline';

    jmodes = 128; % number of modes

    [Thet,Rad] = meshgrid(theta,r);



    %% Initial condition

    r1 = 0; % center of vortex (radius)
    theta1 = 0; % center of vortex (angle)
    circ1 = 1; % vortex circulation

    % r2 = 3*R/8; % center of vortex (radius)
    % theta2 = pi; % center of vortex (angle)
    % circ2 = -1; % vortex circulation

    % r3 = 3*R/8; % center of vortex (radius)
    % theta3 = 0; % center of vortex (angle)
    % circ3 = 1; % vortex circulation


    wf = tanh((R-Rad)/sqrt(2)).*exp(i*7*Thet)... % background
        .*(Rad.*exp(i*circ1*Thet) - r1*exp(i*circ1*theta1))... % 1st vortex
        ./sqrt(Rad.^2 + r1.^2 - 2*r1*Rad.*cos(circ1*(Thet-theta1)) + 1);
    if exist('r2','var')
        wf = wf.*(Rad.*exp(i*circ2*Thet) - r2*exp(i*circ2*theta2))... % 2nd vortex
            ./sqrt(Rad.^2 + r2.^2 - 2*r2*Rad.*cos(circ2*(Thet-theta2)) + 1);
    end
    if exist('r3','var')
        wf = wf.*(Rad.*exp(i*Thet) - r3*exp(i*theta3))... % 3rd vortex
            ./sqrt(Rad.^2 + r3.^2 - 2*r3*Rad.*cos(circ3*(Thet-theta3)) + 1);
    end

    [indi,indj] = find(isinf(wf));
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
    colorbar




    %% Computation of the integration kernel (optional)
    % When the quasi fast Hankel transform is performed, the computation of the
    % so-called "integration kernel" is a long task but does not depend on the
    % input function. It can (has to) be done outside the temporal loop.
    comp_ker = 0
    if comp_ker == 1

        for ii = -nth/2+1:nth/2
            [H,kk,rr,I,KK,RR]=dht([],R,jmodes,ii);
            save(['output/ker_' num2str(floor(ii))],'I','kk','rr','KK','RR')
        end

    end


    t = 0;
    nameout = 'output/psi_00000000';
    save(nameout,'wf','t')
    np = 8; % # of digits following psi_

    % In this file, all the parameters are saved, not just those which will be
    % used later. It can be useful if there is a doubt on the parameters used.
    save('output/misc','r','theta','T','nt','dt','nth','dth','nr','dr','R',...
        'ppskip','Rad','Thet','np','V','comp_ker','interp_sch','wf','t',...
        'split_NL','r1','theta1','circ1')
    if exist('r3','var')
        save('output/misc','r2','theta2','circ2','r3','theta3','circ3','-append')
    elseif exist('r2','var')
        save('output/misc','r2','theta2','circ2','-append')
    end


else % then reload all the parameters and last iteration

    load('output/misc.mat')
    nt = ntn;
    T = dt*nt;
    save('output/misc','r','theta','T','nt','dt','nth','dth','nr','dr','R',...
        'ppskip','Rad','Thet','np','V','comp_ker','interp_sch','wf','t',...
        'split_NL','r1','theta1')
    if exist('r3','var')
        save('output/misc','r2','theta2','circ2','r3','theta3','circ3','-append')
    elseif exist('r2','var')
        save('output/misc','r2','theta2','circ2','-append')
    end

    if size(num2str(nti),2) == np
        nameout = ['output/psi_' num2str(nti)];
    elseif size(num2str(nti),2) == np-1
        nameout = ['output/psi_0' num2str(nti)];
    elseif size(num2str(nti),2) == np-2
        nameout = ['output/psi_00' num2str(nti)];
    elseif size(num2str(nti),2) == np-3
        nameout = ['output/psi_000' num2str(nti)];
    elseif size(num2str(nti),2) == np-4
        nameout = ['output/psi_0000' num2str(nti)];
    elseif size(num2str(nti),2) == np-5
        nameout = ['output/psi_00000' num2str(nti)];
    elseif size(num2str(nti),2) == np-6
        nameout = ['output/psi_000000' num2str(nti)];
    elseif size(num2str(nti),2) == np-7
        nameout = ['output/psi_0000000' num2str(nti)];
    else
        error('Loaded file will be wrong.')
    end
    load(nameout)

end



%% Temporal loop

pp = 1+nti; % indicator of the progression in the loop

tic

deriv_thet = i*(ones(nr,1)*(-nth/2+1:nth/2));

for t = (nti+1)*dt:dt:T

    % 1st step: the linear step, computation of the laplacian part.

    [kt,fr] = obngfft(Thet,wf,2); % fft to isolate every angular mode
    if abs((pp-1)/ppskip - floor((pp-1)/ppskip)) <= 1e-14
        % derivative of wf % theta FOR THE PREVIOUS ITERATION
        [r2,dwfdt] = obngifft(kt,deriv_thet.*fr,2);
    end

    % Hankel transform is easier when performed along r for each mode =>
    % loop.

    for ii = -nth/2+1:nth/2

        load(['output/ker_' num2str(round(ii))])

        % forward Hankel transform
        Fr = interp1(r,fr(:,ii+nth/2),rr,interp_sch,'extrap');
        adht = dht(Fr,RR,KK,I);

        % inverse Hankel transform
        Fr = idht(adht.*exp(-i*0.5*kk.^2*dt),I,KK,RR);
        fr(:,ii+nth/2) = interp1(rr,Fr,r,interp_sch,'extrap');

    end

    [r2,wf] = obngifft(kt,fr,2);
    if abs((pp-1)/ppskip - floor((pp-1)/ppskip)) <= 1e-14
        % the derivatives of wf are added to the previous
        % file which name didnt change since last iteration
        save(nameout,'dwfdt','-append')
    end

    %mass_L = sum(sum(abs(wf).^2.*Rad))*dr*dth

    % 2nd step: NL part

    if split_NL == 0.5 || split_NL == 0 || abs(split_NL-1/3) < 1e-13

        wf = wf.*exp(-i*(V + 0.5*abs(wf).^2)*dt);

    end


    % 3rd step: dissipation

    if split_NL == 0 || abs(split_NL-1/3) < 1e-13

        error('Dissipation not implemented yet.');

    end
    
    %{
    if abs(pp/ppskip - floor(pp/ppskip)) <= 1e-14
        % in this case, an output file is produced
        % name of the outputfile
        if size(num2str(pp),2) == np
            nameout = ['output/psi_' num2str(pp)];
        elseif size(num2str(pp),2) == np-1
            nameout = ['output/psi_0' num2str(pp)];
        elseif size(num2str(pp),2) == np-2
            nameout = ['output/psi_00' num2str(pp)];
        elseif size(num2str(pp),2) == np-3
            nameout = ['output/psi_000' num2str(pp)];
        elseif size(num2str(pp),2) == np-4
            nameout = ['output/psi_0000' num2str(pp)];
        elseif size(num2str(pp),2) == np-5
            nameout = ['output/psi_00000' num2str(pp)];
        elseif size(num2str(pp),2) == np-6
            nameout = ['output/psi_000000' num2str(pp)];
        elseif size(num2str(pp),2) == np-7
            nameout = ['output/psi_0000000' num2str(pp)];
        else
            error('Output names will be wrong.')
        end
        save(nameout,'wf','t')
    end
    %}
    
    pp = pp+1;
    wfi = [wf wf(:,1)];
    %figure()
    colormap(jet(256));
    Z = abs(wfi).^2;
    pcolor(X,Y,Z), shading interp
    %colorbar()
    %plot(X(128,:), Z(128, :));
    %axis equal tight
    xlabel('x')
    ylabel('y')
    colorbar()
    drawnow()
    
end

toc


