% close loop computation for planner

L2 = 1.1556;
mu = 0.012149;
c = zeros(2,1);
c(1) = (1-mu)/(L2 + mu)^3 + mu/(L2 - 1 + mu)^3;
c(2) = 1/(L2 + mu)^3 - 1/(L2 - 1 + mu)^3;
h = zeros(2,1);
h(1) = -1/(2*c(1)*(1-c(1)))*((2-c(1))*4*(L2-2*c(2)*mu*(1-mu))-4*L2);
h(2) = -1/(2*c(1)*(1-c(1)))*(-2*(4*L2-2*c(2)*mu*(1-mu))+4*L2*(1+c(1)));
k = 0.1; 

Omega = 1.8636;
Omega_z = Omega; 
Rotate_w = 0.1;

a = zeros(2,1);
b = zeros(2,1);
a(1) = -1;
a(2) = 1/2;
b(1) = 2;
b(2) = 5/2;
phi = 0;
e   =  0.0549;

x = sym('x',[6 1]);
u = sym('u',[3 1]);
v = sym('v',[3 1]);
z = sym('z', [4 1]);
w = sym('w', [4 1]);

%% linear analytic control model
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

x_1 = x1(1);
x_2 = x1(2);
x_3 = x1(3);
xd_1 = x2(1);
xd_2 = x2(2);

f2 = [x_1 + 2*xd_2 - (mu*(mu + x_1 - 1))/((mu + x_1 - 1)^2 + x_2^2 + x_3^2)^(3/2) + ((mu + x_1)*(mu - 1))/((mu + x_1)^2 + x_2^2 + x_3^2)^(3/2);
    x_2 - 2*xd_1 - (mu*x_2)/((mu + x_1 - 1)^2 + x_2^2 + x_3^2)^(3/2) + (x_2*(mu - 1))/((mu + x_1)^2 + x_2^2 + x_3^2)^(3/2);
    (x_3*(mu - 1))/((mu + x_1)^2 + x_2^2 + x_3^2)^(3/2) - (mu*x_3)/((mu + x_1 - 1)^2 + x_2^2 + x_3^2)^(3/2)];

f2real = -M*x1 - 2*N*x2 - (2*M*x1 + 2*N*x2)*z(1) ...
           - M*x1*z(2) - N*x1*z(3) - ed1*(1 - mu)/absed1^3 ...
           - ed2*mu/absed2^3;

P44 = -2*((x_1+mu)^2+x_2^2+x_3^2)/((x_1+mu)^2+x_2^2+x_3^2)^(5/2) - 2*((x_1-1+mu)^2+x_2^2+x_3^2)/((x_1-1+mu)^2+x_2^2+x_3^2)^(5/2);
P54 = 3*(x_2*(x_1+mu))/((x_1+mu)^2+x_2^2+x_3^2)^(5/2) - 3*(x_2*(x_1-1+mu))/((x_1-1+mu)^2+x_2^2+x_3^2)^(5/2);
P64 = 3*(x_3*(x_1+mu))/((x_1+mu)^2+x_2^2+x_3^2)^(5/2) - 3*x_3*((x_1-1+mu))/((x_1-1+mu)^2+x_2^2+x_3^2)^(5/2);

P = [ 2*x_1 + 2*xd_2, x_1,  x_2, P44;
    2*x_2 - 2*xd_1, x_2, -x_1, P54;
    0,   0,    0, P64];

F = [x2;f2+ P*z + u];

%% Reg feedback
Tz = [b(1)*e 0;0 0;0 -b(1)*e;mu*(1-mu)*a(1)*e 0];
Tr = [k*(1-c(1)+Omega^2)/(2*Omega) 0;0 -k;-k 0];
S2 = [0 -Omega;Omega 0];
Seq = [L2;0;0];
wz = w(1:2);
wr = w(3:4);
pi = [-Tr*wr + Seq;-Tr*S2*wr];

x_1reg = pi(1);
x_2reg = pi(2);
x_3reg = pi(3);
xd_1reg = pi(4);
xd_2reg = pi(5);

f2reg = [x_1reg + 2*xd_2reg - (mu*(mu + x_1reg - 1))/((mu + x_1reg - 1)^2 + x_2reg^2 + x_3reg^2)^(3/2) + ((mu + x_1reg)*(mu - 1))/((mu + x_1reg)^2 + x_2reg^2 + x_3reg^2)^(3/2);
    x_2reg - 2*xd_1reg - (mu*x_2reg)/((mu + x_1reg - 1)^2 + x_2reg^2 + x_3reg^2)^(3/2) + (x_2reg*(mu - 1))/((mu + x_1reg)^2 + x_2reg^2 + x_3reg^2)^(3/2);
    (x_3reg*(mu - 1))/((mu + x_1reg)^2 + x_2reg^2 + x_3reg^2)^(3/2) - (mu*x_3reg)/((mu + x_1reg - 1)^2 + x_2reg^2 + x_3reg^2)^(3/2)];

P44reg = -2*((x_1reg+mu)^2+x_2reg^2+x_3reg^2)/((x_1reg+mu)^2+x_2reg^2+x_3reg^2)^(5/2) - 2*((x_1reg-1+mu)^2+x_2reg^2+x_3reg^2)/((x_1reg-1+mu)^2+x_2reg^2+x_3reg^2)^(5/2);
P54reg = 3*(x_2reg*(x_1reg+mu))/((x_1reg+mu)^2+x_2reg^2+x_3reg^2)^(5/2) - 3*(x_2reg*(x_1reg-1+mu))/((x_1reg-1+mu)^2+x_2reg^2+x_3reg^2)^(5/2);
P64reg = 3*(x_3reg*(x_1reg+mu))/((x_1reg+mu)^2+x_2reg^2+x_3reg^2)^(5/2) - 3*x_3reg*((x_1reg-1+mu))/((x_1reg-1+mu)^2+x_2reg^2+x_3^2)^(5/2);

Preg = [ 2*x_1reg + 2*xd_2reg, x_1reg,  x_2reg, P44reg;
    2*x_2reg - 2*xd_1reg, x_2reg, -x_1reg, P54reg;
    0,   0,    0, P64reg];
%% closed loop

f1cl = x2;
f2cl = f2 + P*z -f2reg - Preg*Tz*wz - Tr*S2^2*wr - v;
f2cl = simplify(f2cl);

f2clreal = f2real - f2reg - Preg*Tz*wz - Tr*S2^2*wr - v;
f2clreal = simplify(f2clreal);
