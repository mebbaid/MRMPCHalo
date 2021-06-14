% simulation of a nominal halo orbit near L2

clc
clear all

% models parameters and init conditions

L2 = 1.1556;
mu = 0.012149;
c = zeros(2,1);
c(1) = (1-mu)/(L2 + mu)^3 + mu/(L2 - 1 + mu)^3;
c(2) = 1/(L2 + mu)^3 - 1/(L2 - 1 + mu)^3;
h = zeros(2,1);
h(1) = -1/(2*c(1)*(1-c(1)))*((2-c(1))*4*(L2-2*c(2)*mu*(1-mu))-4*L2);
h(2) = -1/(2*c(1)*(1-c(1)))*(-2*(4*L2-2*c(2)*mu*(1-mu))+4*L2*(1+c(1)));
k = 0.02; 

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
e   =  0.0549;  % EM rotating system e

% srp parameters
srp = 0; disturbance = 0;
Gsc = 1360.8 ; % approx kw
sped = 300000; % approx m/s
zeta = 0.9252 ;   % prep. on the space-craft


% satellite init position and velocity
u0 = [0;0;0];
% uncomment the following lines if you want to start near the orbit or near
% the moon surface
insertion_error = 0;
% insertion_error = -0.1; 

Seq = [L2 + insertion_error;0;0;0;0;0];   
ho1 = [(-k*(1-c(1)+Omega^2)/(2*Omega))*cos(0);
       k*sin(0);
       k*cos(0)];
   
diffho1 = [(k*(1-c(1)+Omega^2)/2)*sin(Omega*0);
           Omega*k*cos(Omega*0);
           -Omega_z*k*sin(Omega_z*0)];   
       
ho2 = ho1 + [e*h(1)*cos(phi);
             e*h(2)*sin(phi);
             0];
         
diffho2 = diffho1 + [-e*h(1)*sin(phi);
                    e*h(2)*cos(phi);
                    0];       
x0 = Seq + [ho1;diffho1];

% x0 = Seq + [ho2;diffho2];
% x0 = [L2;0;0;0;0;0];
% x0 = [1;0;0;0;0;0]; % starting near Moon surface

%---------------------------------------------------------------------%
%--------------------------- SIMULATION ------------------------------%
%---------------------------------------------------------------------%
rif = 0;
simTime = 200;

out = sim('orbit.slx');
t = 0:10^-3:simTime;

figure 
subplot(2,1,1)
plot3(out.position(:,1), out.position(:,2), out.position(:,3), 'k','LineWidth', 2);
hold on; grid on;
l = legend('Free Evolution in 3D space');
set(l,'Interpreter','Latex'); l.FontSize = 20;
l = xlabel('$x$'); 
set(l,'Interpreter','Latex'); l.FontSize = 20;
l = ylabel('$y$'); 
set(l,'Interpreter','Latex'); l.FontSize = 20;
l = zlabel('$z$'); 
set(l,'Interpreter','Latex'); l.FontSize = 20;

subplot(2,1,2)
plot3(out.ref(:,1), out.ref(:,2), out.ref(:,3), 'r','LineWidth', 2);
hold on; grid on;
l = legend('Reference quasi halo orbit');
set(l,'Interpreter','Latex'); l.FontSize = 20;
l = xlabel('$x$'); 
set(l,'Interpreter','Latex'); l.FontSize = 20;
l = ylabel('$y$'); 
set(l,'Interpreter','Latex'); l.FontSize = 20;
l = zlabel('$z$'); 
set(l,'Interpreter','Latex'); l.FontSize = 20;