function V = V_EX (Freq,EX_char,t)
    T_tot = 1/Freq;%period of a cycle 
    T = T_tot*0.8; % Active piston movement time span
    t_c = t; % mod(abs(t),T_tot); % Cyclic time
    A = EX_char(1,1);
    L = EX_char(1,2);
    V_Dead = EX_char(1,3); 
    V_Start = L*A; 
    V_max = L/(0.5*T);
    if t_c < (0.25*T)
        V = V_Dead;
    elseif (t_c >= 0.25*T)&&(t_c<(0.75*T))
        V = V_Dead + V_max*(t_c-T*0.25)*A;
    elseif (t_c >=(0.75*T))&&(t_c<=(1.25*T))
        V = V_Dead + V_Start - V_max*(t_c-0.75*T)*A;
    else
        V = V_Dead;
    end
end
