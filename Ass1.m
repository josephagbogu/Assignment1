% Name: Agbogu Chidera J.
% Student No: 101056053
% Translation is only considered in the plane, therefore 2 degrees of
% freedom. The following equation holds: %
%
% Where m is 0.26*9.11*10^-31. Effective Mass of the electrons m %
%
%%CONSTANT DECLARATIONS
mass = 0.26*(9.11*10^-31); %% eff mass of Electrons
dt = 10*(10^-15); %sim for 1000 dt's
kb = 1.3806*10^-23;
vth = sqrt((2*300*kb)/(mass)); %thermal velocity, 2D 
nop = 1000; %number of particles
nos = 1000; %number of steps
P1 = 1; %lower bound
P2 = 3; % upper bound
vnet = zeros(1,nop);
x = zeros(1,nop);
y = zeros(1,nop);
vx = zeros(1,nop);
vy = zeros(1,nop);
T = zeros(1,nop);


for t = 1:nop  %initial conditions
    
    x(t) = rand*200*(10^-9);
    y(t) = rand*100*(10^-9);
    vx(t) = (vth)*cos(2*pi*rand);
    vy(t) = (vth)*sin(2*pi*rand);
    vnet(t) = sqrt(vx(t)^2 + vy(t)^2);
end

for s = 1:nos %amount of timesteps 
    
    x(1:nop) = x(1:nop) + (vx(1:nop).*dt); 
    y(1:nop) = y(1:nop) + (vy(1:nop).*dt);

    for n = 1:nop %boundary conditions 
        if (y(n)<=0 || y(n)>=100*10^-9)
           vy(n) = -vy(n);
        end
        if(x(n)<=0)
           x(n) = x(n) + 200*10^-9;
        end
        if(x(n)>=200*10^-9)
          x(n) = x(n) - 200*10^-9;
        end
        
     end

    for f = P1:P2
        colorVec = hsv(5);
        plot(x(f),y(f),'-','color', colorVec(f,:));
    end
    
    T(s) = (mass/(2*kb))*( mean(abs(vx).^2)+ mean(abs(vy).^2) );
    
      axis([0,200*10^-9,0,100*10^-9]); 
      pause(0.00001);
      hold on;
end

figure(2)
h=linspace(0,1000*dt,1000);

plot(h,T)
title('Temp vs. Time');
xlabel('Time (s)'); 
ylabel('Temp (K)'); 
axis([0,10^-14,0,400]);
% 



% COLLISIONS WITH MEAN FREE PATH
% Temperature varies for each scattering of all particles. Average for all
% time steps is rough 300K. More with more particles averge speeds tends to
% vth. Most populous bin in histogram approaches vth. Mean free path and
% mean time between collisions is calculated based on steps/count of % scattering. With more timesteps, MFP and MTBC converges to analytic % result.
clear clf

mass = 0.26*(9.11*10^-31); %effective mass 
dt = 10*(10^-15); %sim for 1000 dt's 
pscat = 1-exp(-(dt)/(0.2*10^(-12))); 
kb = 1.3806*10^-23;
vth = sqrt((2*300*kb)/(mass)); %thermal velocity, 2D
nop = 1000; %no of particles
nos = 1000; %no of steps
P1 = 1; %partition lower bound
P2 = 3; %partition upper bound
count=0;

pcomp = zeros(1,nos);
vnet = zeros(1,nop);
x = zeros(1,nop);
y = zeros(1,nop);
vx = zeros(1,nop);
vy = zeros(1,nop);
T = zeros(1,nop);

for t = 1:nop
   x(t) = rand*200*(10^-9);
   y(t) = rand*100*(10^-9);
   vx(t) = (vth).*randn(1,1)*(1/sqrt(2)); 
   vy(t) = (vth).*randn(1,1)*(1/sqrt(2));
   vnet(t) = sqrt(vx(t)^2 + vy(t)^2);
   %initial conditions
end

for s = 1:nos
    pcomp(s) = rand;
end


for l = 1:nos
    x(1:nop) = x(1:nop) + (vx(1:nop).*dt);
    y(1:nop) = y(1:nop) + (vy(1:nop).*dt);
    
    if(pscat > pcomp(l))
        count = count+1;
        for p = 1:nop
           vx(p) = (vth).*randn(1,1)*(1/sqrt(2)); 
           vy(p) = (vth).*randn(1,1)*(1/sqrt(2));
        end
    end
    
    for n = 1:nop %boundary conditions
        if(y(n)<=0 || y(n)>=100*10^-9)
           vy(n) = -vy(n);
        end
        if(x(n)<=0)
          x(n) = x(n) + 200*10^-9;
        end
        if(x(n)>=200*10^-9)
           x(n) = x(n) - 200*10^-9;
        end
    end
    
    for f=P1:P2
        colorVec = hsv(5);
        plot(x(f),y(f),'-','color', colorVec(f,:));
    end
    
    T(l) = (mass/(2*kb))*( mean(abs(vx).^2)+ mean(abs(vy).^2) );
    
    axis([0,200*10^-9,0,100*10^-9]); 
    pause(0.00001);
    hold on;
end
figure(2); 
h=linspace(0,1000*dt,1000); 
plot(T);
title('Temp vs. Time steps'); 
xlabel('Time Steps'); 
ylabel('Temp (K)');

figure(3);
histogram(vnet,20);
title('Freq vs. Avg Velocity'); 
xlabel('Avg. Velocity (m/s)'); 
ylabel('Freq');

mfp = (nos/count)*dt*mean(vnet); %%Calculation for mean free path
meant = dt*(nos/count); %Mean time between collisions
sprintf('avg. velocity is %0.5e m/s' ,mean(vnet))
sprintf('mean free path is %0.5e m and mean time between collisions is%0.5e s' ,mfp,meant)

