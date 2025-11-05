% Function for based on the adibatic model proposed by :
% I. Urieli and D. M. Berchowitz, Stirling cycle engine analysis. Bristol: A. Hilger, 1984.
% Assumptions:
% System symmetry
% 1 :Ideal gas assumption 
% 2: The pressure is uniform 

function dpdt = AdiabaticDiff(t,Y,Physical_char,Loop_char,WF)

% Variables

p = Y(1); % system pressure
m_c = Y(2); % compressor mass 
W = Y(3); % Work output 
Q_k = Y(4); % Condensor Q_dot  
Q_r = Y(5); % Regenerator Q_dot
Q_h = Y(6);% Heater Q_dot
dpdt = 0*Y; % Answers 


%Physical Characteristics 

Freq = Physical_char(1,1); %System frequency 
M = Physical_char(1,2); %System mass charge 
A_CO = Physical_char(1,3); %Compressor Piston Area
L_CO = Physical_char(1,4); %Compressor Piston Length
V_CO_DEAD = Physical_char(1,5); %Compressor Piston 
V_k = Physical_char(1,6); % Cooler Volume
V_h = Physical_char(1,7); % Heater Volume 
V_r = Physical_char(1,8); % Regen Volume 
%V_p_av = A_CO*L_CO;

CO_char = [A_CO,L_CO,V_CO_DEAD];
EX_char = [A_CO,L_CO,V_CO_DEAD];

% Thermal Loop Characteristics
T_k = Loop_char(1,1); %Heater temperature
T_h = Loop_char(1,2); %Cooler temperature

% Initial values
T_r = (T_h-T_k)/log(T_h/T_k);
Cp = 1006.8;%py.CoolProp.CoolProp.PropsSI('CPMASS', 'T', T_r, 'P', 101325, WF);
Cv = 718.4669;%py.CoolProp.CoolProp.PropsSI('CVMASS', 'T', T_r, 'P', 101325, WF);
gamma = Cp/Cv;
R_1 = 8.314;%py.CoolProp.CoolProp.PropsSI('gas_constant', WF);
Molar = 2.896546000000000e-02;%py.CoolProp.CoolProp.PropsSI('M', WF);
R = R_1/Molar;
DV_c = Vel_CO(Freq,CO_char,t)*A_CO;
DV_e = Vel_EX(Freq,EX_char,t)*A_CO;
V_c  = V_CO(Freq,CO_char,t);
V_e = V_EX(Freq,EX_char,t);
m_k = p*V_k/(R*T_k);
m_r = p*V_r/(R*T_r);
m_h = p*V_h/(R*T_h);
m_e = M-(m_c+m_k+m_h+m_r);
T_c = p*V_c/(R*m_c);
T_e = p*V_e/(R*m_e);


check = Vel_CO(Freq,CO_char,t)<=0;
T_ck = check*T_c + (1-check)*T_k;
%if Vel_CO(Freq,CO_char,t)<=0
  %  T_ck = T_c;
%else
 %   T_ck = T_k;
%end

check = Vel_EX(Freq,EX_char,t)>=0;
T_he = check*T_h + (1-check)*T_e;
% if Vel_EX(Freq,EX_char,t)>=0
%     T_he = T_h;
% else
%     T_he = T_e;
% end

% Equations

Dp = (-gamma*p*(DV_c/T_ck+DV_e/T_he))/(V_c/T_ck+gamma*(V_k/T_k+V_r/T_r +V_h/T_h)+V_e/T_he);
Dm_c = (p*DV_c+V_c*Dp/gamma)/(R*T_ck);
Dm_k = m_k*Dp/p;
Dm_r = m_r*Dp/p;
Dm_h = m_h*Dp/p;
gA_ck = -Dm_c;
gA_kr = gA_ck - Dm_k;
gA_rh = gA_kr - Dm_r;
gA_he = gA_rh - Dm_h;

T_rh = (gA_rh>=0)*T_r + (1-(gA_rh>=0))*T_h;
% if gA_rh >=0
%     T_rh = T_r;
% else
%     T_rh = T_h;
% end
T_kr = (gA_kr>=0)*T_k + (1-(gA_kr>=0))*T_r;
% if gA_kr>=0
%     T_kr = T_k;
% else
%     T_kr = T_r;
% end

DW = p*(DV_c+DV_e);
DQ_k  =  V_k*Dp*Cv/R   - Cp*(T_ck *gA_ck-T_kr*gA_kr);
DQ_r  =  V_r*Dp*Cv/R    - Cp*(T_kr  * gA_kr-T_rh*gA_rh);
DQ_h = (V_h*Dp*Cv/R)  - Cp*(T_rh *gA_rh-T_he*gA_he);

% Return Statements

dpdt(1) = Dp; %Returns pressure differential
dpdt(2) = Dm_c; %Returns compressor mass 
dpdt(3) = DW; % Returns total cycle work
dpdt(4) = DQ_k; %Returns Heat Transfer of the condenser    
dpdt(5) = DQ_r; %Returns Heat Transfer of the regenerator  
dpdt(6) = DQ_h;%Returns Heat Transfer of the heater

end
