%% ADIABATIC STRILING MODEL 
% I. Urieli and D. M. Berchowitz, Stirling cycle engine analysis. Bristol: A. Hilger, 1984.
% Symmetrical assumptions

clear all 
close all 
clc

%Units
Kelvin = 273.15;
R = 287; %For air, for different WF must be changed 
resolution = 1000;
cycles = 3;

% Known Quantities 
WF = 'Air'; % Working fluid 
T_h = 65 + Kelvin; % Heater Temperature
T_k = 10 + Kelvin; % Cooler Temperature 
T_r = (T_h-T_k)/log(T_h/T_k); % Regenerator temperature 
V_k = 8.8e-3; %Cooler Volume 
V_h = 8.8e-3; %Heater Volume 
V_r = 2e-3; % Regenerator Volume, an assumption at this point
Freq = 1.67; % Cycle frequency 
EX_A = 0.00246301; % Expander Piston Area
EX_L = 0.54; %Expnader Swept Length 
%EX_AMAX = 50; %Expander maximum acceleration % Not used 
EX_DEAD = EX_A*0.02;
CO_A =  0.00197213; % Compressor Piston Area % Symmetric assumption
CO_L = EX_A*EX_L/CO_A; % Compressor Length 
CO_DEAD = CO_A*0.02;

%Debug
EX_char = [CO_A,CO_L,CO_DEAD];
CO_char = [CO_A,CO_L,CO_DEAD];

% System charging characterisitcs 
P_charge = 3e5; % System charge pressure 
V_pistons = V_EX(Freq,EX_char,0) + V_CO(Freq,CO_char,0);
M = (V_r+V_h+V_k+V_pistons)*py.CoolProp.CoolProp.PropsSI('D', 'T', T_r, 'P', P_charge, WF);

% Variables needed to solve differential equation 
% system characterisitcs
Physical_char = [Freq,M,CO_A,CO_L,CO_DEAD,V_k,V_h,V_r];
Loop_char = [T_k,T_h];
% Time spans 
T = (1/Freq)*0.8;
epsilon = T*0.25*0.01;
t_1 = linspace(0,(T*0.25)-epsilon,resolution);
t_2 = linspace((T*0.25)+epsilon,(T*0.5)-epsilon,resolution);
t_3 = linspace((T*0.5)+epsilon,(T*0.75)-epsilon,resolution);
t_4 = linspace((T*0.75)+epsilon,T-epsilon,resolution);
t_5 = linspace((T*1.0)+epsilon,(1.25*T)-epsilon,resolution); 
% Initial conditions 
p_init = P_charge;
m_c_init = CO_A*CO_L*py.CoolProp.CoolProp.PropsSI('D', 'T', T_r, 'P', P_charge, WF);
W_init = 0; % Independent variable
Q_k_init = 0; % Independent variable
Q_r_init = 0; % Independent variable
Q_h_init = 0; % Independent variable
init_1 = [p_init,m_c_init,W_init,Q_k_init,Q_r_init,Q_h_init];

%% Differential Equations 

[t1,Y1] = ode113(@(t1,Y1)AdiabaticDiff(t1,Y1,Physical_char,Loop_char,WF),t_1,init_1);
init_2 = Y1(end,:);
[t2,Y2] = ode113(@(t2,Y2)AdiabaticDiff(t2,Y2,Physical_char,Loop_char,WF),t_2,init_2);
init_3 = Y2(end,:);
[t3,Y3] = ode113(@(t3,Y3)AdiabaticDiff(t3,Y3,Physical_char,Loop_char,WF),t_3,init_3);
init_4 = Y3(end,:);
[t4,Y4] = ode113(@(t4,Y4)AdiabaticDiff(t4,Y4,Physical_char,Loop_char,WF),t_4,init_4);
init_5 = Y4(end,:);
[t5,Y5] = ode113(@(t5,Y5)AdiabaticDiff(t5,Y5,Physical_char,Loop_char,WF),t_5,init_5);

%% Processing the data

% Mass
M_1_2 = trapz (t1,Y1(:,3));
M_2_3 = trapz(t2,Y2(:,3));
M_3_4 = trapz(t3,Y3(:,3));
M_4_5 = trapz(t4,Y4(:,3));
M_5_1 = trapz(t5,Y5(:,3));


% Work
W_1_2 = trapz (t1,Y1(:,3));
W_2_3 = trapz(t2,Y2(:,3));
W_3_4 = trapz(t3,Y3(:,3));
W_4_5 = trapz(t4,Y4(:,3));
W_5_1 = trapz(t5,Y5(:,3));
% W_1_2 = trapz (t1,Y1(:,3));
% W_2_3 = trapz(t2,Y2(:,3));
% W_3_4 = trapz(t3,Y3(:,3));
% W_4_5 = trapz(t4,Y4(:,3));
% W_5_1 = trapz(t5,Y5(:,3));
W_net = W_1_2+W_2_3+W_3_4+W_4_5+W_5_1;

% Heater  
Q_1_2 = trapz (t1,Y1(:,6));
Q_2_3 = trapz(t2,Y2(:,6));
Q_3_4 = trapz(t3,Y3(:,6));
Q_4_5 = trapz(t4,Y4(:,6));
Q_5_1 = trapz(t5,Y5(:,6));
Qh_net = Q_1_2+Q_2_3+Q_3_4+Q_4_5+Q_5_1;

% Cooler 
Qk_1_2 = trapz (t1,Y1(:,4));
Qk_2_3 = trapz(t2,Y2(:,4));
Qk_3_4 = trapz(t3,Y3(:,4));
Qk_4_5 = trapz(t4,Y4(:,4));
Qk_5_1 = trapz(t5,Y5(:,4));
Qk_net = Qk_1_2+Qk_2_3+Qk_3_4+Qk_4_5+Qk_5_1;

%% Data Plotting 

% Pressure plot
plot(t1,Y1(:,1),t2,Y2(:,1),t3,Y3(:,1),t4,Y4(:,1),t5,Y5(:,1))
ylabel('Pressure [pa]')
xlabel('time [s]')
title('Pressure per time in the engine')

% Mass in compressor plot 
figure 
plot(t1,Y1(:,2),t2,Y2(:,2),t3,Y3(:,2),t4,Y4(:,2),t5,Y5(:,2))
ylabel('Mass [kg]')
xlabel('time [s]')
title('Mass in compressor piston per time')

% T-S diagrams 
S_sytem = zeros(5,resolution);
Densities = zeros(5,resolution);
for i = 1:resolution
    Desnity = Y1(i,2)/V_CO(Freq,CO_char,t1(i,1));
    S_sytem(1,i) = py.CoolProp.CoolProp.PropsSI('S', 'D',Desnity , 'P', Y1(i,1), WF);
end 
for i = 1:resolution
    Desnity = Y2(i,2)/V_CO(Freq,CO_char,t2(i,1));
    S_sytem(2,i) = py.CoolProp.CoolProp.PropsSI('S', 'D',Desnity , 'P', Y2(i,1), WF);
end 
for i = 1:resolution
    Desnity = Y3(i,2)/V_CO(Freq,CO_char,t3(i,1));
    S_sytem(3,i) = py.CoolProp.CoolProp.PropsSI('S', 'D',Desnity , 'P', Y3(i,1), WF);
end 
for i = 1:resolution
    Desnity = Y4(i,2)/V_CO(Freq,CO_char,t4(i,1));
    S_sytem(4,i) = py.CoolProp.CoolProp.PropsSI('S', 'D',Desnity , 'P', Y4(i,1), WF);
end 
for i = 1:resolution
    Desnity = Y5(i,2)/V_CO(Freq,CO_char,t5(i,1));
    S_sytem(5,i) = py.CoolProp.CoolProp.PropsSI('S', 'D',Desnity , 'P', Y5(i,1), WF);
end 
figure
plot(t1,S_sytem(1,:),t2,S_sytem(2,:),t3,S_sytem(3,:),t4,S_sytem(4,:),t5,S_sytem(5,:));

% Q_h and Q_k and W
figure 
plot(t1,Y1(:,6),t2,Y2(:,6),t3,Y3(:,6),t4,Y4(:,6),t5,Y5(:,6))
hold on
plot(t1,Y1(:,4),t2,Y2(:,4),t3,Y3(:,4),t4,Y4(:,4),t5,Y5(:,4))
hold on
plot(t1,Y1(:,3),t2,Y2(:,3),t3,Y3(:,3),t4,Y4(:,3),t5,Y5(:,3))
hold on 
%plot(t1,Y1(:,5),t2,Y2(:,5),t3,Y3(:,5),t4,Y4(:,5),t5,Y5(:,5))
