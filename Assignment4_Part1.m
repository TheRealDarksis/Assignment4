%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        ELEC 4700 A - Winter 2022                        %
%                              Assignment 4                               %
%        Author: Julie-Anne Chaine             Student #: 101104568       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
ny = 30; nx = 20;
W = 20e-9; L = 30e-9;
V = zeros(nx, ny);
vMap = zeros(nx,ny);

V_range = linspace(0.1,10,10);          

TB = 0.7e-7;
BB = 0.3e-7;
RB = 1.2e-7;
LB = 0.8e-7;

% Constants:
mo = 9.1093837015e-31;              % Electron rest mass in kg
mn = 0.26*mo;                       % Effective mass of electrons
k = 1.38e-23;                       % Boltzmann's constant 
W = 200e-9;                         % Nominal region width in m
L = 100e-9;                         % Nominal region length in m
tau_mn = 0.2e-12;                   % Mean time between collisions
q = 1.602e-19;                      % Electron charge
concentration = 1e19; % in m

%Variables
T = 300;                            % Room temperature in K
elecpop = 1000;                     % Number of particles to simulate
dt = -2.5e-14;                         % Time step in s
x = zeros(elecpop,1);               % Electron y positions
y = zeros(elecpop,1);               % Electron x positions
oldx = zeros(elecpop,1);            % Previous electron x position
oldy = zeros(elecpop,1);            % Previous electron y position
Vth = sqrt(2*k*T/mn);               % Thermal velocity in m/s 
samplepop = 10;                     % # of particles to plot
samp = randi(elecpop,samplepop,1);  % Random particles to observe
iter = 1000;                         % # of iterations 
Pscatter = 1 - exp(-dt/tau_mn);     % Scattering probability
ax = zeros(elecpop,1);
ay = zeros(elecpop,1);

% Initializing positions and velocities 
x(:,1) = rand(elecpop,1)*W;
y(:,1) = rand(elecpop,1)*L;
vx = Vth*randn(elecpop,1)/sqrt(2);
vy = Vth*randn(elecpop,1)/sqrt(2); 
    for z = 1:elecpop   % Generating points on the left side
       x(z,1) = 0.2e-7 + (0.8e-7 - 0.2e-7)*rand();
    end
    for z = elecpop/2+1:elecpop     % Generating points on the right
       x(z,1) = 1.2e-7 + (2e-7 - 1.2e-7)*rand();        
    end 

    for j = 1:length(V_range) 
        cMap = ones(nx, ny);
        for i = 8:12
            for indexy = 1:7
                cMap(i,indexy) = 0.5;
            end
            for indexy = 22:30
                cMap(i,indexy) = 0.5;
            end
        end

        for i = 1:nx 
            for jj = 1:ny
                nxm = jj + (i-2)*ny;
                nxp = jj + i*ny;
                nyp = jj + 1 + (i-1)*ny;
                nym = jj - 1 + (i-1)*ny;
                n = jj + (i-1)*ny;
                if i == 1           % Left
                    G(n, n) = cMap(i,jj);
                    F(n) = V_range(j);
                elseif i == nx      % Right
                    G(n, n) = cMap(i,jj);
                    F(n) = 0;
                elseif jj == 1       % Bottom     
                    G(n, nxm) = cMap(i-1,jj);         
                    G(n, nxp) = cMap(i+1,jj);          
                    G(n, nyp) = cMap(i,jj+1); 
                    G(n, n) = -(G(n, nxp)+G(n, nxm)+G(n, nyp));
                elseif jj == ny      % Top
                    G(n, nxm) = cMap(i-1,jj);         
                    G(n, nxp) = cMap(i+1,jj);                   
                    G(n, nym) = cMap(i,jj-1); 
                    G(n, n) = -(G(n, nxp)+G(n, nxm)+G(n, nym));
                else
                    G(n, nxm) = cMap(i-1,jj);         
                    G(n, nxp) = cMap(i+1,jj);          
                    G(n, nyp) = cMap(i,jj+1);         
                    G(n, nym) = cMap(i,jj-1);  
                    G(n, n) = -(G(n, nxp)+G(n, nxm)+G(n, nym)+G(n, nyp));
                end
            end
        end
        V = G\F';

        for i = 1:nx
            for indexy = 1:ny
                n = indexy + (i-1)*ny;
                vMap(i,indexy) = V(n);
            end
        end
        [Ey,Ex] = gradient(vMap);
        Ex = -Ex*nx/W;
        Ey = -Ey*ny/L;
        
        for index = 1:iter
            Ybin = discretize(y,ny);
            Xbin = discretize(x,nx);
            for b = 1:elecpop
                Ex_particle = Ex(Xbin(b), Ybin(b));
                Ey_particle = Ey(Xbin(b), Ybin(b));
                Fx = Ex_particle*q;
                ax(b) = Fx/mn;
                Fy = Ey_particle*q;
                ay(b) = Fy/mn;
            end
            
            %Creating box reflection, 0.4e-7 bottle neck
            oldMBox = oldx>LB & oldx<RB & oldy<TB & oldy>BB;
            MBox = x>LB & x<RB & (y>TB | y<BB);
            vy(oldMBox & MBox) = -vy(oldMBox & MBox);
            vx(MBox & ~oldMBox) = -vx(MBox & ~oldMBox);
            
            oldx = x;
            oldy = y;

            x = oldx + vx*dt; % Previous x position + delta L + acc
            y = oldy + vy*dt; % Previous y position + delta L + acc
            vx = vx + ax*dt;
            vy = vy + ay*dt;


            PartScatter = Pscatter > rand(elecpop,1);
            vx(PartScatter) = Vth*randn(sum(PartScatter),1)/sqrt(2);
            vy(PartScatter) = Vth*randn(sum(PartScatter),1)/sqrt(2);

            %Creating box reflection, 0.6e-7 bottle neck
            oldMBox = oldx>LB & oldx<RB & oldy<TB & oldy>BB;
            MBoxb = x>LB & x<RB & y<BB;
            MBoxt = x>LB & x<RB & y>TB;
            MBox = MBoxb | MBoxt;
            vy(oldMBox & MBox) = -vy(oldMBox & MBox);
            vx(MBox & ~oldMBox) = -vx(MBox & ~oldMBox);

            % If the particle enters the box:
            x(oldx>RB & MBox & ~oldMBox) = RB;
            x(oldx<LB & MBox & ~oldMBox) = LB;
            y(oldMBox & MBoxb) = BB;
            y(oldMBox & MBoxt) = TB;

            oldx(x<0) = W; % Making the particles on the left boundary, appear on the right
            oldx(x>W) = 0; % Making the particles on the right boundary, appear on the left
            x(x<0) = x(x<0) + W; % All points passed the left get pushed to the other side
            x(x>W) = x(x>W) - W; % All points passed the right get pushed to the other side
            vy(y>L) = -vy(y>L); % Particle Y direction gets flipped if it hits the top
            vy(y<0) = -vy(y<0); % Particle Y direction gets flipped
         
        end
        
        V_calc = sqrt( vx.^2 + vy.^2 );
        Vmean = mean(V_calc);
        Vmean_x = mean(vx);
        J(j) = mean(q*concentration*Vmean_x*L); % In a 2D cross-section we use L
    end
plot(V_range,J)
title('Current vs Voltage behaviour')
xlabel('Voltage (V)')
ylabel('Current (A)')
Resistors = polyfit(J, V_range,1)
