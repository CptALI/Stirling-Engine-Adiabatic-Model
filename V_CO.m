function V = V_CO(Freq,CO_char,t)
    T_tot = 1/Freq;%period of a cycle 
    T = T_tot*0.8; % Active piston movement time span
    t_c = t; % mod(abs(t),T_tot); % Cyclic time
    A = CO_char(1,1);
    L = CO_char(1,2);
    V_Dead = CO_char(1,3);
    V_Start = L*A; 
    V_max = L/(0.5*T);
    if t_c < (0.5*T)
        V = V_Start - V_max*t_c*A + V_Dead;
    elseif (t_c>=0.5*T)&&(t_c<(0.75*T))
        V = V_Dead;
    elseif (t_c>=(0.75*T))&&(t_c<(1.25*T))
        V = V_Dead + V_max*(t_c-0.75*T)*A;
    else
        V = V_Dead;
    end
end