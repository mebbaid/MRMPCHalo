clc
clear all

%% set case to be simulated
delta = 0.1; % adjust for hours (0.15 is one hour) (0.01 = 4 minutes)
emulation  = delta; 
% delta_b = delta/2;
saturation = 1; % set to one to incorporate saturation on the control.
sat_constraint = 1; % set to one to include saturation as as a constraint in MPC formulation
disturbance = 1;  % set to one to incoporate disturbances
srp         = 1; % 
delay = 0; % put to one to include effect of delay
if saturation == 1
    satValue = 0.55; % adjust this value to saturate inputs
else
    satValue = inf;
end


%% 
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
% Omega_z = Omega; 
Omega_z = 1.79; %out-of-plane frequency
Rotate_w = 0.1;

a = zeros(2,1);
b = zeros(2,1);
a(1) = -1;
a(2) = 1/2;
b(1) = 2;
b(2) = 5/2;
phi = 0;
e   =  0.0549;  % EM rotating system e
% e   = 0;

% srp parameters
Gsc = 1360.8 ; % approx kw
sped = 300000; % approx m/s
zeta = 0.9252 ;   % prep. on the space-craft


% satellite init position and velocity
u0 = [0;0;0];
x0 = [L2;0;0;0;0;0];
% x0 = [1;0;0;0;0;0]; % starting near Moon surface
%uncomment the following lines if you want to start near the orbit or near
%the moon surface
% insertion_error = 0;
insertion_error = -0.1; 
Seq = [L2 + insertion_error;0;0;0;0;0];   

ho1 = [(-k*(1-c(1)+Omega^2)/(2*Omega))*cos(0);
       k*sin(0);
       k*cos(0)];
   
diffho1 = [(k*(1-c(1)+Omega^2)/2)*sin(Omega*0);
           Omega*k*cos(Omega*0);
           -Omega_z*k*sin(Omega_z*0)];
      
% x0 = Seq + [ho1;diffho1];

% feedback linearization and pole placement if activated
pole = [-2 -1.5 -3 -2.5 -3.1 -2.6];
A21 = [1+2*c(1) 0 0; 0 1-c(1) 0; 0 0 0-c(1)];
A22 = [0 2 0; -2 0 0; 0 0 0];
% A = [zeros(3) eye(3); zeros(3) zeros(3)];
A = [zeros(3) eye(3); A21 A22];
G = [zeros(3); eye(3)];
K = place(A,G,pole);

%%
simTime = 20; % set to 4380 for long term 6 months station keeping
ref_select = 1;  % set to 1 for L2 orbit, set to 2 for to consider also effects of eccentricity
t = 0:10^-3:simTime;
ts = 0:delta:simTime;

out1 = sim('orbit.slx');

Omega_z = Omega; 

out2 = sim('orbit.slx');



figure('Name','Reference_Type')
subplot(1,2,1);
plot3(out1.ref(:,1), out1.ref(:,2), out1.ref(:,3), 'b', 'LineWidth', 2)
hold on; grid on;
scatter3(L2,0,0,'b','diamond');
l = xlabel('$x$');
set(l,'Interpreter','Latex');
l.FontSize = 15;
l = ylabel('$y$');
set(l,'Interpreter','Latex');
l.FontSize = 15;
l = zlabel('$z$');
set(l,'Interpreter','Latex');
l.FontSize = 15;
l = legend('quasi Halo orbit', '$L_2$ point');
set(l,'Interpreter','Latex');
l.FontSize = 20;


subplot(1,2,2);
plot3(out2.ref(:,1), out2.ref(:,2), out2.ref(:,3), 'b', 'LineWidth', 2)
hold on; grid on;
scatter3(L2,0,0,'b','diamond');
l = xlabel('$x$');
set(l,'Interpreter','Latex');
l.FontSize = 15;
l = ylabel('$y$');
set(l,'Interpreter','Latex');
l.FontSize = 15;
l = zlabel('$z$');
set(l,'Interpreter','Latex');
l.FontSize = 15;
l = legend('quasi Halo orbit', '$L_2$ point');
set(l,'Interpreter','Latex');
l.FontSize = 20;