clc
clear all

%% set case to be simulated
delta = 0.1; % adjust for hours (0.15 is one hour) (0.01 = 4 minutes)
% delta_b = delta/2;
saturation = 0; % set to one to incorporate saturation on the control.
sat_constraint = 0; % set to one to include saturation as as a constraint in MPC formulation
disturbance = 1;  % set to one to incoporate disturbances
srp         = 1; % solar radiation pressure
emulation =delta; % put emulation = delta to simulate emulated control for FL
delay = 0; % put to one to include effect of delay
if saturation == 1
    satValue = 15; % adjust this value to saturate inputs
else
    satValue = inf;
end

Ts = delta;
timescale = 6.5; % scaling factor such that each s in simulation is an hour
distanceScale = 384400; % distance between two primaries
errorScale = distanceScale/100;  % error in km
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
Gsc = 1.3608 ; % approx kw
sped = 300000; % approx m/s
zeta = 0.9252 ;   % prep. on the space-craft

% satellite init position and velocity
u0 = [0;0;0];
% x0 = [L2;0;0;0;0;0];
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
      
x0 = Seq + [ho1;diffho1];

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
reffl = ans.xr;
n_reffl = ans.n_ref;
yfl = ans.y*distanceScale;
zfl = ans.z;
vfl = ans.v;
ufl = ans.u;

xreg = ans.xreg;
refreg = ans.refreg;
yreg = ans.yreg*distanceScale;
zreg = ans.zreg;
ureg = ans.ureg;

efl = ans.error*errorScale;
e_rms_fl = ans.e_rms*errorScale;
deltaVfl = ans.DeltaV;



figure('Name','Feedback linearization');

subplot(3,2,1);
l = title('x-z plot');
set(l,'Interpreter','Latex');
plot(reffl(:,1)*distanceScale, reffl(:,3)*distanceScale, 'k', 'LineWidth', 1.5);
hold on; grid on;
plot(yfl(:,1), yfl(:,3), 'r', 'LineWidth', 1.5);
l = legend('$x_r$-$z_r$ Nominal reference', '$x$-$z$ FL actual trajectory');
set(l,'Interpreter','Latex');
l = xlabel('y-z plot (km)'); 
l.FontSize = 18;


subplot(3,2,2)
l = title('y-z plot');
set(l,'Interpreter','Latex');
plot(reffl(:,2)*distanceScale, reffl(:,3)*distanceScale, 'k', 'LineWidth', 1.5);
hold on; grid on;
plot(yfl(:,2), yfl(:,3), 'r', 'LineWidth', 1.5);
l = legend('$y_r$-$z_r$ Nominal reference', '$y$-$z$ FL actual trajectory');
set(l,'Interpreter','Latex');
l = xlabel('y-z plot (km)'); 
l.FontSize = 18;

subplot(3,2,3);
l = title('x-y plot');
set(l,'Interpreter','Latex');
plot(reffl(:,1)*distanceScale, -reffl(:,2)*distanceScale, 'k', 'LineWidth', 1.5);
hold on; grid on;
plot(yfl(:,1), -yfl(:,2), 'r', 'LineWidth', 1.5);
l = legend('$x_r$-$y_r$ Nominal reference', '$x$-$-y$ FL actual trajectory');
set(l,'Interpreter','Latex');
l = xlabel('x-y plot (km)'); 
l.FontSize = 18;

subplot(3,2,4);
l = title('x-z plot');
set(l,'Interpreter','Latex');
plot(reffl(:,1)*distanceScale, reffl(:,3)*distanceScale, 'k', 'LineWidth', 1.5);
hold on; grid on;
plot(yfl(:,1), yfl(:,3), 'r', 'LineWidth', 1.5);
l = legend('$x_r$-$z_r$ Nominal reference', '$x$-$z$ FL actual trajectory ');
set(l,'Interpreter','Latex');
l = xlabel('x-z plot (km)'); 
l.FontSize = 18;

subplot(3,2,5);
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


subplot(3,2,6);
l = title('Primer disturbances');
set(l,'Interpreter','Latex');
plot(t*timescale, zfl(:,1), 'k', 'LineWidth', 1.5);
hold on; grid on;
plot(t*timescale, zfl(:,2), 'r', 'LineWidth', 1.5);
plot(t*timescale, zfl(:,3), 'b', 'LineWidth', 1.5);
plot(t*timescale, zfl(:,4), 'g', 'LineWidth', 1.5);
l = legend('$z_1(t)$', '$z_2(t)$', '$z_3(t)$','$z_4(t)$');
set(l,'Interpreter','Latex');
l = xlabel('Time (h)'); 
l.FontSize = 18;





figure('Name','Nonlinear regulation');

subplot(3,2,1);
l = title('x-z plot');
set(l,'Interpreter','Latex');
plot(refreg(:,1)*distanceScale, refreg(:,3)*distanceScale, 'k', 'LineWidth', 1.5);
hold on; grid on;
plot(yreg(:,1), yreg(:,3), 'r', 'LineWidth', 1.5);
l = legend('$x_r$-$z_r$ Nominal reference', '$x$-$z$ FL actual trajectory');
set(l,'Interpreter','Latex');
l = xlabel('y-z plot (km)'); 
l.FontSize = 18;


subplot(3,2,2)
l = title('y-z plot');
set(l,'Interpreter','Latex');
plot(refreg(:,2)*distanceScale, refreg(:,3)*distanceScale, 'k', 'LineWidth', 1.5);
hold on; grid on;
plot(yreg(:,2), yreg(:,3), 'r', 'LineWidth', 1.5);
l = legend('$y_r$-$z_r$ Nominal reference', '$y$-$z$ FL actual trajectory');
set(l,'Interpreter','Latex');
l = xlabel('y-z plot (km)'); 
l.FontSize = 18;

subplot(3,2,3);
l = title('x-y plot');
set(l,'Interpreter','Latex');
plot(refreg(:,1)*distanceScale, -refreg(:,2)*distanceScale, 'k', 'LineWidth', 1.5);
hold on; grid on;
plot(yreg(:,1), -yreg(:,2), 'r', 'LineWidth', 1.5);
l = legend('$x_r$-$y_r$ Nominal reference', '$x$-$-y$ FL actual trajectory');
set(l,'Interpreter','Latex');
l = xlabel('x-y plot (km)'); 
l.FontSize = 18;

subplot(3,2,4);
l = title('x-z plot');
set(l,'Interpreter','Latex');
plot(refreg(:,1)*distanceScale, refreg(:,3)*distanceScale, 'k', 'LineWidth', 1.5);
hold on; grid on;
plot(yreg(:,1), yreg(:,3), 'r', 'LineWidth', 1.5);
l = legend('$x_r$-$z_r$ Nominal reference', '$x$-$z$ FL actual trajectory ');
set(l,'Interpreter','Latex');
l = xlabel('x-z plot (km)'); 
l.FontSize = 18;

subplot(3,2,5);
l = title('Controls');
set(l,'Interpreter','Latex');
plot(t*timescale, ureg(:,1), 'k', 'LineWidth', 1.5);
hold on; grid on;
plot(t*timescale, ureg(:,2), 'r', 'LineWidth', 1.5);
plot(t*timescale, ureg(:,3), 'b', 'LineWidth', 1.5);
l = legend('$u_1(t)$', '$u_2(t)$', '$u_3(t)$');
set(l,'Interpreter','Latex');
l = xlabel('Time (h)'); 
l.FontSize = 18;


subplot(3,2,6);
l = title('Primer disturbances');
set(l,'Interpreter','Latex');
plot(t*timescale, zreg(:,1), 'k', 'LineWidth', 1.5);
hold on; grid on;
plot(t*timescale, zreg(:,2), 'r', 'LineWidth', 1.5);
plot(t*timescale, zreg(:,3), 'b', 'LineWidth', 1.5);
plot(t*timescale, zreg(:,4), 'g', 'LineWidth', 1.5);
l = legend('$z_1(t)$', '$z_2(t)$', '$z_3(t)$','$z_4(t)$');
set(l,'Interpreter','Latex');
l = xlabel('Time (h)'); 
l.FontSize = 18;
