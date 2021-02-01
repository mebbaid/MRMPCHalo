function xr = ref_gen(t,k,c,Omega, Omega_z, Seq, Ts)
%REF_GEN Summary of this function goes here
%   Detailed explanation goes here
    ho1 = [(-k*(1-c(1)+Omega^2)/(2*Omega))*cos(t);
           k*sin(t);
           k*cos(t)];

    diffho1 = [(k*(1-c(1)+Omega^2)/2)*sin(Omega*t);
               Omega*k*cos(Omega*t);
               -Omega_z*k*sin(Omega_z*t)];
    
    xr1 = Seq + [ho1;diffho1];       
           
    ho11 = [(-k*(1-c(1)+Omega^2)/(2*Omega))*cos(Omega*(t+Ts));
           k*sin(Omega*(t+Ts));
           k*cos(Omega_z*(t+Ts))];

    diffho11 = [(k*(1-c(1)+Omega^2)/2)*sin(Omega*(t+Ts));
               Omega*k*cos(Omega*(t+Ts));
               -Omega_z*k*sin(Omega_z*(t+Ts))];

    xr2 = Seq + [ho11;diffho11];

    xr = [xr1 xr2];       
           
end

