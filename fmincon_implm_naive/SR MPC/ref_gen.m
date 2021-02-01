function xr = ref_gen(xr,np,k,c,Omega,Omega_z,delta,Seq,t)
%REF_GEN Summary of this function goes here
%   Detailed explanation goes here
% for i = 1:np  
%         ho1 = [(-k*(1-c(1)+Omega^2)/(2*Omega))*cos(i*delta);
%                 k*sin(i*delta);
%                 k*cos(i*delta)];
% 
%          diffho1 = [(k*(1-c(1)+Omega^2)/2)*sin(Omega*i*delta);
%                Omega*k*cos(Omega*i*delta);
%                -Omega_z*k*sin(Omega_z*i*delta)];
% 
%         xr(:,i) = Seq + [ho1;diffho1];
% end
        ho1 = [(-k*(1-c(1)+Omega^2)/(2*Omega))*cos(t);
                k*sin(t);
                k*cos(t)];

         diffho1 = [(k*(1-c(1)+Omega^2)/2)*sin(Omega*t);
               Omega*k*cos(Omega*t);
               -Omega_z*k*sin(Omega_z*t)];

        xr = Seq + [ho1;diffho1];
end

