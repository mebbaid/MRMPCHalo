clc
clear all

%% set case to be simulated
delta = 0.05; % adjust for hours (0.15 is one hour) (0.01 = 4 minutes)
% delta_b = delta/2;
saturation = 0; % set to one to incorporate saturation on the control.
sat_constraint = 0; % set to one to include saturation as as a constraint in MPC formulation
disturbance = 0;  % set to one to incoporate disturbances
srp         = 0; % solar radiation pressure
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
% insertion_error = -0.3; 
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
ts = 0:delta:simTime;

% MR inversion
sim('HaloSim_SD.slx'); %uncomment to check MR inversion
xmr = ans.x;
refmr = ans.xr;
n_refmr = ans.n_ref;
ymr = ans.y*distanceScale;
zmr = ans.z;
vmr = ans.v;
umr = ans.u;
emr = ans.error*errorScale;
e_rms_mr = ans.e_rms;
deltaVmr = ans.DeltaV;

veloIncr_mr = sum(deltaVmr(round(simTime*timescale/2)))/21;


figure('Name','MR SD regulation');

subplot(2,2,1);
l = title('x-z plot');
set(l,'Interpreter','Latex');
plot3(refmr(:,1)*distanceScale, refmr(:,2)*distanceScale,refmr(:,3)*distanceScale, 'k', 'LineWidth', 1.5);
hold on; grid on;
plot3(ymr(:,1), ymr(:,2),ymr(:,3), 'r', 'LineWidth', 1.5);
scatter3(L2*distanceScale,0,0,'b','diamond');
l = legend('$x(t), y(t), z(t)$ Nominal reference', '$x(t), y(t), z(t)$ MR actual trajectory', 'L2 point');
set(l,'Interpreter','Latex');
l = xlabel('x-z plot (km)'); 
l.FontSize = 18;


subplot(2,2,2)
l = title('Controls');
set(l,'Interpreter','Latex');
plot(t*timescale, umr(:,1), 'k', 'LineWidth', 1.5);
hold on; grid on;
plot(t*timescale, umr(:,2), 'r', 'LineWidth', 1.5);
plot(t*timescale, umr(:,3), 'b', 'LineWidth', 1.5);
l = legend('$u_1(t)$', '$u_2(t)$', '$u_3(t)$');
set(l,'Interpreter','Latex');
l = xlabel('Time (h)'); 
l.FontSize = 18;


subplot(2,2,3);
l = title('Norm of the error');
set(l,'Interpreter','Latex');
plot(ts*timescale, e_rms_mr, 'r', 'LineWidth', 1.5);
hold on; grid on;
l = legend('MR SD regulation $\|e(t)\|$');
set(l,'Interpreter','Latex');
l = xlabel('Time (h)'); 
l.FontSize = 18;
hold off;

subplot(2,2,4);
l = title('Control effort magnitutde');
set(l,'Interpreter','Latex');
plot(t*timescale, deltaVmr, 'k', 'LineWidth', 1.5);
hold on; grid on;
l = legend('MR SD regulation  $\|u\|(t)$');
set(l,'Interpreter','Latex');
l = xlabel('Time (h)'); 
l.FontSize = 18;
hold off;
