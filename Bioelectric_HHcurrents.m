%graphing currents of AP 
%Defining constants from variables given in lecture 
gNa = 120;
gK = 36;
gL = 0.3;

Na_out = 490;
K_out = 20;
Na_in = 50;
K_in= 400;

Nernst_leak = -50;
Cm = 1;
k =1;

%initial values 
%Vr 
y0(1,1) = -62; 
%h
y0(2,1) = 0.55; 
%m
y0(3,1) = 0.059; 
%n
y0(4,1) = 0.33; 

%From equation 
Ek = -70.51; 
Ena = 53.72; 

dt = [0,40]; %time of integration
options = odeset('RelTol', 1e-4, 'AbsTol', [1e-8 1e-8 1e-8 1e-8],'MaxStep', 0.01);

%Integrating model with ode45
[t,y] = ode45('hh_diff_eq', dt,y0, options);

Ik = gK.*(((y(:,4).^4).*((y(:,1))-Ek)));
Ina = gNa.*((y(:,3).^3).*(y(:,2)).*(y(:,1)-Ena));

%attempting to plot variables 
figure
subplot(3,1,1);
plot(t,y(:,1));
xlabel('Time (ms)');
ylabel('Voltage (mV)');
title('Voltage Time Series');

subplot(3,1,2);
plot(t,y(:,2));
hold all
xlabel('Time(ms)');
ylabel('Probability of Opening');
plot(t,y(:,3));
xlabel('Time (ms)');
ylabel('Probability of Opening');
plot(t,y(:,4));
xlabel('Time (ms)');
ylabel('Probability of Opening');
title('Gating Time Series');
legend('h gate','m gate','n gate');

subplot(3,1,3);
plot(t,Ina);
hold all
xlabel('Time (ms)');
ylabel('Current (mA)');
title('Current vs. Time');
plot(t, Ik);
legend('Ina','Ik');