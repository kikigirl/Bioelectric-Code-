%using function command 

function Bioelectric_hw2(ext_time,dist,AA) 

%%%dist= distance in cm, AA = amplitude of input, ext_time = runtime of potential  

%%%including given parameters 

%%%converting to kilohms 
          %units 
Re = .35; %kohm-cm 
Ri = .11; %kohm-cm 
Cm = 2.5; %uF/cm^2
a_radius = 10^-3; %cm
diameter = 2*a_radius; %cm
width = 10^-4; %cm
delta_x = 0.1; %cm
g_Na = 1445; %mS/cm2
g_L = 128; %mS/cm2
E_na = 115; %mV
E_l = -0.01; %mV

%%%new parameters 
dt = 0.001; %time interval in ms
time = ext_time;

%%%instantiating time vector 
t = 0:dt:time;
t_first = length(t);

%%%instantiating voltage vectors 
V=zeros(91,t_first);
Ve=zeros(91,t_first);
Vrest= zeros(91,t_first);
dVdt=zeros(91,t_first); %for change in membrane potential

%%%instantiating variable changes 
m=zeros(91,t_first);
h=zeros(91,t_first);
Im=zeros(91,t_first);  %for membrane current 
dmdt=zeros(91,t_first);
dhdt=zeros(91,t_first);

%%%using given parameter equations 
alpha_m0=97/(1+exp(31/5.3));   %%%simplified for easier computations
beta_m0=alpha_m0/exp(-23.8/4.17);
beta_h0=15.6/(1+exp(24/10));
alpha_h0=beta_h0/exp(-5.5/5);

%%%including parameter vectors 
m(:,1)=alpha_m0/(alpha_m0+beta_m0);
h(:,1)=alpha_h0/(alpha_h0+beta_h0);

%%%first loop to fill the Ist (stimulus) current 
for t = 1:t_first-1
    
    if t < 150 
        Ist = AA;
    else 
        Ist = 0;
    end 
    
    %%%second loop to update Ve (membrane potential) 
    for x=1:91
        Ve(x,t)=Re*Ist/(4*pi*sqrt((0.1*x-4.5)^2+dist^2)); %referenced from lecture notes 
    end
    
    %%%third loop to update outputs of equations 
    for x=2:90 
        Vm=V(x,t)-Vrest(x,t);  %measuring change in membrane potential
        alpha_m=(97+0.363*Vm)/(1+exp((31-Vm)/5.3));
        beta_m=alpha_m/exp((Vm-23.8)/4.17);
        beta_h=15.6/(1+exp((24-Vm)/10));
        alpha_h=beta_h/exp((Vm-5.5)/5);
        dmdt(x,t)=-1*(alpha_m+beta_m)*m(x,t)+alpha_m;
        dhdt(x,t)=-1*(alpha_h+beta_h)*h(x,t)+alpha_h;
        Im(x,t)=g_Na*m(x,t)^2*h(x,t)*(Vm-E_na)+g_L*(Vm-E_l);
        dVdt(x,t)=(1/Cm)*(-Im(x,t)+(diameter*delta_x)/(4*Ri*width)*(((V(x-1,t)-2*V(x,t)+V(x+1,t))/(delta_x^2))+((Ve(x-1,t)-2*Ve(x,t)+Ve(x+1,t))/(delta_x^2))));
    end 
    
    %%%using Euler's method for approximations 
    if t~t_first;
    V(:,t+1)=V(:,t)+dVdt(:,t)*dt;
    m(:,t+1)=m(:,t)+dmdt(:,t)*dt;
    h(:,t+1)=h(:,t)+dhdt(:,t)*dt;
    end
end

%%%plotting surface potential of Vm(x,t) 
x=1:91;
t=0:dt:ext_time;

surf(t,x,V,'edgecolor', 'none')

title(['Bi-Directional Propagating Action Potential from ',num2str(dist),'cm away'])
xlabel('Time (ms)')
ylabel('Axon Position (mm)')
zlabel('Vm (mV)')

%%%plotting video of bi-directional AP 
figure 
vid=VideoWriter('APBi_direct_video','MPEG-4');
open(vid);
for i= 1:20:t_first
   plot(x,V(x,i))
  title(['Movie of a Bi-Directional Propagating Action Potential from ',num2str(dist),'cm away'])
    xlabel('Axon Position(mm)')
    ylabel('Vm(mV)')
    axis([0,91,-20,100])
    mov=getframe;
    writeVideo(vid,mov)
end
close(vid)
end