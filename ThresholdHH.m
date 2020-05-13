function [dy] = hh_diff_eq(t,y)
%Unload passed variables
V = y(1,1);
h = y(2,1);
m = y(3,1);
n = y(4,1);

%for #2 square pulse duration of .25ms
% if t > 5 && t < 5.25
% I = 20;
% else
% I = 0;
% end 

%for problem #3a
% if t > 5 && t < 5.25
%     I = 20; % need increased Is to generate 2nd AP
% elseif t > 23 && t < 23.25 %18ms before 2nd AP depolarized
%     I = 20;
% else
%     I = 0;
% end

%for problem 3b
% if t > 5 && t < 5.25
%     I = 20*3;  %tripling the magnitude 
% elseif t > 16 && t < 16.25 
%     I = 20*3;
% else
%     I = 0;
% end

%for problem 4
if t > 5 && t < 5.25
    I = -20*3;
else
    I = 0;
end

gNa= 120;
gK = 36;
gL = 0.3;

Na_out= 490;
K_out = 20;
Na_in = 50;
K_in= 400;

Nernst_leak = -50;
Cm = 1;
k = 1;

Ek = -70.51; 
Ena = 53.72; 

Vdiff = V+62; %-62 = Rest Potential

a_m = ((.1*(25-Vdiff))/(exp((25-Vdiff)/10)-1));
a_h = .07*exp(-Vdiff/20);
a_n = (.01*(10-Vdiff))/(exp((10-Vdiff)/10)-1);

b_n = .125*exp(-Vdiff/80);
b_m = 4*exp(-Vdiff/18);
b_h = 1/(exp((30-Vdiff)/10)+1);

dy(1,1) = (1/Cm)*(I-(gNa*(m^3)*h*(V-Ena))-(gK*(n^4)*(V-Ek))-(gL*(V-Nernst_leak))); dy(2,1) = ((-(a_h+b_h)*h) +a_h)*k;
dy(3,1) = ((-(a_m+b_m)*m) +a_m)*k;
dy(4,1) = ((-(a_n+b_n)*n) +a_n)*k;

end