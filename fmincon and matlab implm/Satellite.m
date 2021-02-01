function dx = Satellite(x,u,Ts)
%SATELLITE Summary of this function goes here
%   Detailed explanation goes here

L2 = 1.1556;
mu = 0.012149;
c = zeros(2,1);
c(1) = (1-mu)/(L2 + mu)^3 + mu/(L2 - 1 + mu)^3;
c(2) = 1/(L2 + mu)^3 - 1/(L2 - 1 + mu)^3;
h = zeros(2,1);
h(1) = -1/(2*c(1)*(1-c(1)))*((2-c(1))*4*(L2-2*c(2)*mu*(1-mu))-4*L2);
h(2) = -1/(2*c(1)*(1-c(1)))*(-2*(4*L2-2*c(2)*mu*(1-mu))+4*L2*(1+c(1)));
Omega = 1.8636;
Omega_z = Omega; 

z = zeros(4,1);

x1 = x(1:3);
x2 = x(4:6);
M = [-1 0 0;
    0 -1 0;
    0 0 0];

N = [0 -1 0;
    1 0 0;
    0 0 0];

d1r = [-mu - z(4)/(1-mu);0;0];
d2r = [1-mu + z(4)/mu;0;0];

ed1 = x1 - d1r;
absed1 = sqrt(ed1(1)^2+ed1(2)^2+ed1(3)^2);
ed2 = x1 - d2r;
absed2 = sqrt(ed2(1)^2+ed2(2)^2+ed2(3)^2);

f1 = x2;
f2 = -M*x1 - 2*N*x2 - (2*M*x1 + 2*N*x2)*z(1) ...
           - M*x1*z(2) - N*x1*z(3) - ed1*(1 - mu)/absed1^3 ...
           - ed2*mu/absed2^3;
F = [f1; f2];
G = [zeros(3); eye(3)];
dx = F + G*u;

end

