% main functions gahering different simulations and cases 

clc
clear all

%% set case to be simulated
delta = 0.05; % adjust for hours (0.15 is one hour) (0.01 = 4 minutes)
% delta_b = delta/2;
saturation = 0; % set to one to incorporate saturation on the control.
sat_constraint = 0; % set to one to include saturation as as a constraint in MPC formulation
disturbance = 0;  % set to one to incoporate disturbances
srp         = 0; % 
emulation = delta; % put emulation = delta to simulate emulated control for FL
delay = 0; % put to one to include effect of delay
if saturation == 1
    satValue = 15; % adjust this value to saturate inputs
else
    satValue = inf;
end

r = 0.1; % control penalty if consistent penalty on three controls is required

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

%% mpc parameters
nx = 6;
ny = 6;  % outputs are x,y directly
nu = 3;
np = 4;
nc = 4;
Ts = delta;
timescale = 6.5; % scaling factor such that each s in simulation is an hour
distanceScale = 384400; % distance between two primaries
errorScale = distanceScale/1000;  % error in km
accScale = 1000; % accelarations in m/s^2
options = optimoptions('fmincon','Algorithm','sqp');
% Nonlinear MPC controller
nlobj = nlmpc(nx, ny, nu);
nlobj.Ts = Ts;
nlobj.PredictionHorizon = np;
nlobj.ControlHorizon = nc; 
nlobj.Model.StateFcn = "Satellite";
nlobj.Model.IsContinuousTime = true;
nlobj.Model.OutputFcn = @(x,u,Ts) [x(1);x(2);x(3); x(4);x(5);x(6)];
nlobj.Model.NumberOfParameters = 1;
nlobj.Weights.OutputVariables = [10 10 10 1 1 1];
% nlobj.Weights.OutputVariables = [1 1 1 1 1 1];
nlobj.Weights.ManipulatedVariablesRate = r*[1 1 1];% try to play with weights 
if sat_constraint == 1
    nlobj.ManipulatedVariables(1).Max = satValue;
    nlobj.ManipulatedVariables(2).Max = satValue;
    nlobj.ManipulatedVariables(3).Max = satValue;
    nlobj.ManipulatedVariables(1).Min = -satValue;
    nlobj.ManipulatedVariables(2).Min = -satValue;
    nlobj.ManipulatedVariables(3).Min = -satValue;
end
validateFcns(nlobj,x0,u0,[],{Ts}); 
createParameterBus(nlobj,'MPCHalo/Nonlinear MPC Controller','myBusObject',{Ts});

%% Simulation 
simTime = 10; % set to 4380 for long term 6 months station keeping
ref_select = 1;  % set to 1 for L2 orbit, set to 2 for to consider also effects of eccentricity
t = 0:10^-3:simTime;

sim('MPCHalo.slx');
xmpc = ans.x;
refmpc = ans.xr;
n_refmpc = ans.n_ref;
ympc = ans.y*distanceScale;
ysr = ans.ysr*distanceScale;
zmpc = ans.z;
vmpc = ans.v;
umpc = ans.u;
e_rms_mpc = ans.e_rms;
norm_umpc = ans.DeltaV;
mrmpcstatus = ans.mrmpcstatus;
% veloIncr_mrmpc = sum(deltaVmpc(round(simTime*timescale/2)))/21;




figure('Name','MR MPC');

subplot(2,2,1);
l = title('Proposed MR MPC');
set(l,'Interpreter','Latex');
plot3(ympc(:,1),ympc(:,2), ympc(:,3), 'r', 'LineWidth', 3);
hold on; grid on;
% plot3(refmpc(:,1)*distanceScale,refmpc(:,2)*distanceScale, refmpc(:,3)*distanceScale, 'b--', 'LineWidth', 1.5);
plot3(n_refmpc(:,1)*distanceScale,n_refmpc(:,2)*distanceScale, n_refmpc(:,3)*distanceScale, 'k', 'LineWidth', 3);
scatter3(L2*distanceScale,0,0,'b','diamond');
% plot3(xe1, xe2, xe3, 'k', 'LineWidth', 2);
l = xlabel('$x$');
set(l,'Interpreter','Latex');
l = ylabel('$y$');
set(l,'Interpreter','Latex');
l = zlabel('$z$');
set(l,'Interpreter','Latex');
l = legend('$x(t), y(t), z(t)$- MR MPC trajectory','$x(t), y(t), z(t)$- nominal reference', 'L2 point' );
set(l,'Interpreter','Latex');
set(l,'Interpreter','Latex');


subplot(2,2,2);
l = title('Controls');
set(l,'Interpreter','Latex');
plot(t*timescale, umpc(:,1), 'k', 'LineWidth', 1.5);
hold on; grid on;
plot(t*timescale, umpc(:,2), 'r', 'LineWidth', 1.5);
plot(t*timescale, umpc(:,3), 'b', 'LineWidth', 1.5);
l = legend('MR MPC $u_1(t)$', 'MR MPC $u_2(t)$', 'MR MPC $u_3(t)$');
set(l,'Interpreter','Latex');
l = xlabel('Time (h)'); 
l.FontSize = 18;


subplot(2,2,3);
l = title('Norm of the error');
set(l,'Interpreter','Latex');
plot(t*timescale, e_rms_mpc, 'r', 'LineWidth', 1.5);
hold on; grid on;
l = legend('MR MPC $\|e_{RMS}(t)\|$ km');
set(l,'Interpreter','Latex');
l = xlabel('Time (h)'); 
l.FontSize = 18;
hold off;

subplot(2,2,4);
l = title('Control effort magnitutde');
set(l,'Interpreter','Latex');
plot(t*timescale, norm_umpc, 'k', 'LineWidth', 1.5);
hold on; grid on;
l = legend('MR MPC  $\|u\|(t)$');
set(l,'Interpreter','Latex');
l = xlabel('Time (h)'); 
l.FontSize = 18;
hold off;

