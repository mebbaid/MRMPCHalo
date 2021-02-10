clc
clear all

%% set case to be simulated
delta = 0.05; % adjust for hours (0.15 is one hour) (0.01 = 4 minutes)
% delta_b = delta/2;
saturation = 0; % set to one to incorporate saturation on the control.
disturbance = 0;  % set to one to incoporate disturbances
srp         = 0; % solar radiation pressure
emulation = delta; % put emulation = delta to simulate emulated control for FL
delay = 0; % put to one to include effect of delay
if saturation == 1
    satValue = 15; % adjust this value to saturate inputs
else
    satValue = inf;
end

Ts = delta;
timescale = 6.5; % scaling factor such that each s in simulation is an hour
distanceScale = 384400; % distance between two primaries
errorScale = distanceScale/1000;  % error in km
accScale = 1000; % accelarations in m/s^2

%% models parameters and init conditions
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
% insertion_error = -0.1; 
% Seq = [L2 + insertion_error;0;0;0;0;0];   
% 
% ho1 = [(-k*(1-c(1)+Omega^2)/(2*Omega))*cos(0);
%        k*sin(0);
%        k*cos(0)];
%    
% diffho1 = [(k*(1-c(1)+Omega^2)/2)*sin(Omega*0);
%            Omega*k*cos(Omega*0);
%            -Omega_z*k*sin(Omega_z*0)];
%       
% x0 = Seq + [ho1;diffho1];

% feedback linearization and pole placement if activated
pole = [-2 -1.5 -3 -2.5 -3.1 -2.6];
A21 = [1+2*c(1) 0 0; 0 1-c(1) 0; 0 0 0-c(1)];
A22 = [0 2 0; -2 0 0; 0 0 0];
% A = [zeros(3) eye(3); zeros(3) zeros(3)];
A = [zeros(3) eye(3); A21 A22];
G = [zeros(3); eye(3)];
K = place(A,G,pole);

simTime = 15; % set to 4380 for long term 6 months station keeping
ref_select = 1;  % set to 1 for L2 orbit, set to 2 for to consider also effects of eccentricity
t = 0:10^-3:simTime;


% feedback linearization
sim('HaloSim.slx');    
% plots
xfl = ans.x;
reffl = ans.xr*distanceScale;
yfl = ans.y*distanceScale;
zfl = ans.z;
vfl = ans.v;
ufl = ans.u;

xreg = ans.xreg;
refreg = ans.refreg*distanceScale;
yreg = ans.yreg*distanceScale;
zreg = ans.zreg;
ureg = ans.ureg;

e_rms_fl = ans.e_rms;
e_rms_reg = ans.e_rms_reg;

norm_ufl = ans.DeltaV;
norm_ureg = ans.DeltaV_reg;


figure('Name','Feedback linearization');

subplot(2,2,1);
l = title('3D plot of trajectory under FL');
set(l,'Interpreter','Latex');
plot3(reffl(:,1), reffl(:,2), reffl(:,3), 'k', 'LineWidth', 1.5);
hold on; grid on;
plot3(yfl(:,1),yfl(:,2), yfl(:,3), 'r', 'LineWidth', 1.5);
l = legend('$x(t), y(t), z(t)$- Nominal reference','$x(t), y(t), z(t)$- Feedback linearization trajectory', 'L2 point' );
set(l,'Interpreter','Latex');
l = xlabel('y-z plot (km)'); 
l.FontSize = 18;


subplot(2,2,2)
l = title('Controls');
set(l,'Interpreter','Latex');
plot(t*timescale, ufl(:,1), 'k', 'LineWidth', 1.5);
hold on; grid on;
plot(t*timescale, ufl(:,2), 'r', 'LineWidth', 1.5);
plot(t*timescale, ufl(:,3), 'b', 'LineWidth', 1.5);
l = legend('$u_1(t)$', '$u_2(t)$', '$u_3(t)$');
set(l,'Interpreter','Latex');
l = xlabel('Time (h)'); 
l.FontSize = 18;

subplot(2,2,3);
l = title('Norm of the error');
set(l,'Interpreter','Latex');
plot(t*timescale, e_rms_fl, 'r', 'LineWidth', 1.5);
hold on; grid on;
l = legend('Feedback linearization $\|e(t)\|$ km');
set(l,'Interpreter','Latex');
l = xlabel('Time (h)'); 
l.FontSize = 18;
hold off;


subplot(2,2,4);
l = title('Control effort magnitutde');
set(l,'Interpreter','Latex');
plot(t*timescale, norm_ufl, 'k', 'LineWidth', 1.5);
hold on; grid on;
l = legend('Feedback linearization $\|u\|(t)$');
set(l,'Interpreter','Latex');
l = xlabel('Time (h)'); 
l.FontSize = 18;
hold off;


figure('Name','Nonlinear regulation');

subplot(2,2,1);
l = title('3D plot of trajectory under regulation');
set(l,'Interpreter','Latex');
plot3(refreg(:,1), refreg(:,2), refreg(:,3), 'k', 'LineWidth', 1.5);
hold on; grid on;
plot3(yreg(:,1),yreg(:,2), yreg(:,3), 'r', 'LineWidth', 1.5);
l = legend('$x(t), y(t), z(t)$- Nominal reference','$x(t), y(t), z(t)$- Nonlinear regulaiton trajectory', 'L2 point' );
set(l,'Interpreter','Latex');
l = xlabel('y-z plot (km)'); 
l.FontSize = 18;

subplot(2,2,2)
l = title('Controls');
set(l,'Interpreter','Latex');
plot(t*timescale, ureg(:,1), 'k', 'LineWidth', 1.5);
hold on; grid on;
plot(t*timescale, ureg(:,2), 'r', 'LineWidth', 1.5);
plot(t*timescale, ureg(:,3), 'b', 'LineWidth', 1.5);
l = legend('Nonlinear regulation $u_1(t)$', 'Nonlinear regulation $u_2(t)$', 'Nonlinear regulation $u_3(t)$');
set(l,'Interpreter','Latex');
l = xlabel('Time (h)'); 
l.FontSize = 18;

subplot(2,2,3);
l = title('Norm of the error');
set(l,'Interpreter','Latex');
plot(t*timescale, e_rms_reg, 'r', 'LineWidth', 1.5);
hold on; grid on;
l = legend('Nonlinear regulation $\|e(t)\|$ km');
set(l,'Interpreter','Latex');
l = xlabel('Time (h)'); 
l.FontSize = 18;
hold off;


subplot(2,2,4);
l = title('Control effort magnitutde');
set(l,'Interpreter','Latex');
plot(t*timescale, norm_ureg, 'k', 'LineWidth', 1.5);
hold on; grid on;
l = legend('Nonlinear regulation $\|u\|(t)$');
set(l,'Interpreter','Latex');
l = xlabel('Time (h)'); 
l.FontSize = 18;
hold off;