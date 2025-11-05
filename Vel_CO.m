% Step function velocity of the compressore piston 

function Vel = Vel_CO(Freq,CO_char,t)
    T_tot = 1/Freq;%period of a cycle 
    T = T_tot*0.8; % Active piston movement time span
    t_c = mod(abs(t),T_tot); % Cyclic time
    A = CO_char(1,1);
    L = CO_char(1,2);
    V_max = L/(0.5*T);
    if t_c < (0.5*T)
        Vel = - V_max;
    elseif ((0.5*T)<t_c)&&(t_c<(0.75*T))
        Vel = 0;
    elseif (t_c>(0.75*T))&&(t_c<(1.25*T))
        Vel =V_max;
    else
        Vel = 0;
    end

end