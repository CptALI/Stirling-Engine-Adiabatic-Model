% Function for expander piston 

function Vel = Vel_EX(Freq,EX_char,t)
        T_tot = 1/Freq;%period of a cycle 
        T = T_tot*0.8; % Active piston movement time span
        %t_c = mod(abs(t),T_tot); % Cyclic time
        t_c = t;
        A = EX_char(1,1);
        L = EX_char(1,2);
        V_max = L/(0.5*T);
        if t_c < (0.25*T)
            Vel = 0;
        elseif (t_c>=0.25*T)&&(t_c<(0.75*T))
            Vel = V_max;
        elseif (t_c>=(0.75*T))&&(t_c<=(1.25*T))
            Vel =-V_max;
        else
            Vel = 0;
        end
end
