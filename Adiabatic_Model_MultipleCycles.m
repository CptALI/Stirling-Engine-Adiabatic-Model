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
cycles = 50;

% Known Quantities 
WF = 'Air'; % Working fluid 
T_h = 265 + Kelvin; % Heater Temperature
T_k = 10 + Kelvin; % Cooler Temperature 
T_r = (T_h-T_k)/log(T_h/T_k); % Regenerator temperature 
V_k = 8.8e-5; %Cooler Volume 
V_h = 8.8e-5; %Heater Volume 
V_r = 2e-5; % Regenerator Volume, an assumption at this point
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
t_array = zeros(cycles*5,resolution);
for i = 1:cycles  
    cycle_time = (i-1)*(1/Freq);
    t_array((i-1)*5+1,:) = linspace(0+cycle_time,cycle_time+(T*0.25)-epsilon,resolution);
    t_array((i-1)*5+2,:)  = linspace((T*0.25)+epsilon+cycle_time,cycle_time+(T*0.5)-epsilon,resolution);
    t_array((i-1)*5+3,:)  = linspace((T*0.5)+epsilon+cycle_time,cycle_time+(T*0.75)-epsilon,resolution);
    t_array((i-1)*5+4,:)  = linspace((T*0.75)+epsilon+cycle_time,cycle_time+T-epsilon,resolution);
    t_array((i-1)*5+5,:)  = linspace((T*1.0)+epsilon+cycle_time,cycle_time+(1.25*T)-epsilon,resolution);    
end

% Initial conditions first cycle
p_init = P_charge;
m_c_init = CO_A*CO_L*py.CoolProp.CoolProp.PropsSI('D', 'T', T_r, 'P', P_charge, WF);
W_init = 0; % Independent variable
Q_k_init = 0; % Independent variable
Q_r_init = 0; % Independent variable
Q_h_init = 0; % Independent variable
inits = zeros(cycles*5,6);
inits(1,:) = [p_init,m_c_init,W_init,Q_k_init,Q_r_init,Q_h_init];

% Result Cell Arrays  
Y_result = cell(cycles,5); 
t_result = cell(cycles,5);

%% Differential Equations 
% using a nested loop to cycle through 5 states and different cycles 

for i = 1:cycles % cycles loop 
    for j = 1:5 % Processes loop
       [t_result{i,j},Y_result{i,j}] = ode113(@(t1,Y1)AdiabaticDiff(t1,Y1,Physical_char,Loop_char,WF), ...
           t_array((i-1)*5+j,:),inits((i-1)*5+j,:));
       if j<5
           inits((i-1)*5+j+1,:) = Y_result{i,j}(end,:);
       end
    end
    if i<cycles
        inits(i*5+1,:) = Y_result{i,5}(end,:);
    end
end

%% Data Plotting 

% Pressure plot
pressure = zeros(1,resolution*cycles*5);
times = zeros(1,resolution*cycles*5);
for i = 1:cycles
    for j =1:5
        pressure(1,((i-1)*5*resolution+(j-1)*resolution+1):((i-1)*5*resolution+j*resolution)) = Y_result{i,j}(:,1);
        times(1,((i-1)*5*resolution+(j-1)*resolution+1):((i-1)*5*resolution+j*resolution)) = t_result{i,j}(:,1);
    end
end

plot(times,pressure)
ylabel('Pressure [pa]')
xlabel('time [s]')
title('Pressure per time in the Engine')

% Mass in compressor plot 
masses = zeros(1,resolution*cycles*5);
for i = 1:cycles
    for j =1:5
        masses(1,((i-1)*5*resolution+(j-1)*resolution+1):((i-1)*5*resolution+j*resolution)) = Y_result{i,j}(:,2);
    end
end
figure 
plot(times,masses)
ylabel('Mass [kg]')
xlabel('time [s]')
title('Mass in compressor piston per time')

% work plot
works = zeros(1,resolution*cycles*5);
times = zeros(1,resolution*cycles*5);
for i = 1:cycles
    for j =1:5
        works(1,((i-1)*5*resolution+(j-1)*resolution+1):((i-1)*5*resolution+j*resolution)) = Y_result{i,j}(:,3);
        times(1,((i-1)*5*resolution+(j-1)*resolution+1):((i-1)*5*resolution+j*resolution)) = t_result{i,j}(:,1);
    end
end
figure
plot(times,works)
ylabel('Work [J]')
xlabel('time [s]')
title('Work per time in the Engine')

% work plot
Q_regen = zeros(1,resolution*cycles*5);
times = zeros(1,resolution*cycles*5);
for i = 1:cycles
    for j =1:5
        Q_regen(1,((i-1)*5*resolution+(j-1)*resolution+1):((i-1)*5*resolution+j*resolution)) = Y_result{i,j}(:,5);
        times(1,((i-1)*5*resolution+(j-1)*resolution+1):((i-1)*5*resolution+j*resolution)) = t_result{i,j}(:,1);
    end
end
figure
plot(times,Q_regen)
ylabel('Q_regen [J]')
xlabel('time [s]')
title('Q_{regen} vs time in the Engine')

% work plot
Q_heater = zeros(1,resolution*cycles*5);
times = zeros(1,resolution*cycles*5);
for i = 1:cycles
    for j =1:5
        Q_heater(1,((i-1)*5*resolution+(j-1)*resolution+1):((i-1)*5*resolution+j*resolution)) = Y_result{i,j}(:,6);
        times(1,((i-1)*5*resolution+(j-1)*resolution+1):((i-1)*5*resolution+j*resolution)) = t_result{i,j}(:,1);
    end
end
figure
plot(times,Q_heater)
ylabel('Q_{heater} [J]')
xlabel('time [s]')
title('Q_{heater} vs time in the Engine')