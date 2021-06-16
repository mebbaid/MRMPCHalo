% main functions gahering different simulations and cases 

clc
clear all
addpath('./MR MPC simulink implm');
addpath('./SR MPC simulink implm');
addpath('C:\Users\mohd-\Desktop\MRMPCHalo\feedback linearization and regulation');

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

planner_type = 0; % select planner model.
r = 0.0; % control penalty if consistent penalty on three controls is required

%% models parameters and init conditions
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
% Omega_z = 1.79; %out-of-plane frequency
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

%% mpc parameters
nx = 6;
ny = 6;  % outputs are x,y directly
nu = 3;
np = 6;
nc = np;
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
% nlobj.Weights.OutputVariables = [1 1 1 0 0 0;10 10 10 1 1 1];  % scenario 2 
nlobj.Weights.OutputVariables = [1 1 1 1 1 1];
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



% SR MPC Block settings

nlobjsr = nlmpc(nx, ny, nu);
nlobjsr.Ts = Ts;
nlobjsr.PredictionHorizon = np;
nlobjsr.ControlHorizon = nc; 
nlobjsr.Model.StateFcn = "Satellite";
nlobjsr.Model.IsContinuousTime = true;
nlobjsr.Model.OutputFcn = @(x,u,Ts) [x(1);x(2);x(3); x(4);x(5);x(6)];
nlobjsr.Model.NumberOfParameters = 1;
% nlobj.Weights.OutputVariables = [1 1 1 0 0 0;10 10 10 1 1 1];  % scenario 2 
nlobjsr.Weights.OutputVariables = [1 1 1 1 1 1];
nlobjsr.Weights.ManipulatedVariablesRate = r*[1 1 1];% try to play with weights 
if sat_constraint == 1
    nlobjsr.ManipulatedVariables(1).Max = satValue;
    nlobjsr.ManipulatedVariables(2).Max = satValue;
    nlobjsr.ManipulatedVariables(3).Max = satValue;
    nlobjsr.ManipulatedVariables(1).Min = -satValue;
    nlobjsr.ManipulatedVariables(2).Min = -satValue;
    nlobjsr.ManipulatedVariables(3).Min = -satValue;
end
validateFcns(nlobjsr,x0,u0,[],{Ts}); 
createParameterBus(nlobjsr,'MPCHalo/Nonlinear MPC Controller','myBusObject',{Ts});




%% Simulation 
simTime = 25; % set to 4380 for long term 6 months station keeping
ref_select = 1;  % set to 1 for L2 orbit, set to 2 for to consider also effects of eccentricity
t = 0:10^-3:simTime;
ts = 0:delta:simTime;

sim('MPCHalo.slx');

xmpc = ans.x;
refmpc = ans.xr;
n_refmpc = ans.n_ref;
% ympc = ans.y*distanceScale;
% ysr = ans.ysr*distanceScale;
ympc = ans.y;
ysr = ans.ysr;
zmpc = ans.z;
vmpc = ans.v;
umpc = ans.u;
e_rms_mpc = ans.e_rms;
norm_umpc = ans.DeltaV;
mrmpcstatus = ans.mrmpcstatus;
% veloIncr_mrmpc = sum(deltaVmpc(round(simTime*timescale/2)))/21;
clear ans;

sim('HaloMPCMatlab.slx');
xsr = ans.x;
refsr = ans.xr;
% ysr = ans.y*distanceScale;
ysr = ans.y;
zsr = ans.z;
vsr = ans.v;
usr= ans.u;
esr = ans.error;
e_rms_sr = ans.e_rms;
deltaVsr = ans.DeltaV;
srmpcstatus = ans.srmpcstatus;
veloIncr_sr = sum(deltaVsr(round(simTime*timescale/2)))/21;

clear ans;


% feedback linearization
sim('HaloSim.slx');    
% plots
xfl = ans.x;
% reffl = ans.xr*distanceScale;
% yfl = ans.y*distanceScale;
reffl = ans.xr;
yfl = ans.y;
zfl = ans.z;
vfl = ans.v;
ufl = ans.u;

xreg = ans.xreg;
% refreg = ans.refreg*distanceScale;
% yreg = ans.yreg*distanceScale;
refreg = ans.refreg;
yreg = ans.yreg;
zreg = ans.zreg;
ureg = ans.ureg;

e_rms_fl = ans.e_rms;
e_rms_reg = ans.e_rms_reg;

e_fl = ans.error2;
e_reg = ans.error1;

norm_ufl = ans.DeltaV;
norm_ureg = ans.DeltaV_reg;

clear ans;


figure('Name','Comparasions');
subplot(3,3,1);
plot3(reffl(:,1), reffl(:,2), reffl(:,3), 'k', 'LineWidth', 1.5);
hold on; grid on;
mr = plot3(ympc(:,1),ympc(:,2), ympc(:,3), 'b', 'LineWidth', 1);
sr = plot3(ysr(:,1),ysr(:,2), ysr(:,3), 'r','LineWidth', 1);
% reg = plot3(yreg(:,1),yreg(:,2), yreg(:,3), 'color','#D95319', 'LineWidth', 1);
% fl = plot3(yfl(:,1),yfl(:,2), yfl(:,3), 'g', 'LineWidth', 1);
scatter3(L2,0,0,'b','diamond');
xlim([1.1 1.22]);
ylim([-0.1 0.1]);
zlim([-0.2 0.12]);
l = xlabel('$x$');
set(l,'Interpreter','Latex');
l.FontSize = 15;
l = ylabel('$y$');
set(l,'Interpreter','Latex');
l.FontSize = 15;
l = zlabel('$z$');
set(l,'Interpreter','Latex');
l.FontSize = 15;

subplot(3,3,2);
plot3(reffl(:,1), reffl(:,2), reffl(:,3), 'k', 'LineWidth', 1.5);
hold on; grid on;
mr = plot3(ympc(:,1),ympc(:,2), ympc(:,3), 'b', 'LineWidth', 1);
% sr = plot3(ysr(:,1),ysr(:,2), ysr(:,3), 'r','LineWidth', 1);
reg = plot3(yreg(:,1),yreg(:,2), yreg(:,3), 'color','#D95319', 'LineWidth', 1);
% fl = plot3(yfl(:,1),yfl(:,2), yfl(:,3), 'g', 'LineWidth', 1);
scatter3(L2,0,0,'b','diamond');
xlim([1.1 1.22]);
ylim([-0.1 0.1]);
zlim([-0.2 0.12]);
l = xlabel('$x$');
set(l,'Interpreter','Latex');
l.FontSize = 15;
l = ylabel('$y$');
set(l,'Interpreter','Latex');
l.FontSize = 15;
l = zlabel('$z$');
set(l,'Interpreter','Latex');
l.FontSize = 15;


subplot(3,3,3);
plot3(reffl(:,1), reffl(:,2), reffl(:,3), 'k', 'LineWidth', 1.5);
hold on; grid on;
mr = plot3(ympc(:,1),ympc(:,2), ympc(:,3), 'b', 'LineWidth', 1);
% sr = plot3(ysr(:,1),ysr(:,2), ysr(:,3), 'r','LineWidth', 1);
% reg = plot3(yreg(:,1),yreg(:,2), yreg(:,3), 'color','#D95319', 'LineWidth', 1);
fl = plot3(yfl(:,1),yfl(:,2), yfl(:,3), 'g', 'LineWidth', 1);
scatter3(L2,0,0,'b','diamond');
xlim([1.1 1.22]);
ylim([-0.1 0.1]);
zlim([-0.2 0.12]);
l = xlabel('$x$');
set(l,'Interpreter','Latex');
l.FontSize = 15;
l = ylabel('$y$');
set(l,'Interpreter','Latex');
l.FontSize = 15;
l = zlabel('$z$');
set(l,'Interpreter','Latex');
l.FontSize = 15;

subplot(3,3,[4,6]);
plot(ts*timescale, e_rms_mpc, 'b', 'LineWidth', 1.5);
hold on; grid on;
plot(t*timescale, e_rms_sr, 'r', 'LineWidth', 1.5);
plot(t*timescale, e_rms_reg,  'color','#D95319','LineWidth', 1.5);
plot(t*timescale, e_rms_fl, 'g','LineWidth', 1.5);
ylim([-0.1 0.2]);
l = xlabel('$Time (h)$');
set(l,'Interpreter','Latex');
l.FontSize = 15;
l = ylabel('Tracking error');
set(l,'Interpreter','Latex');
l.FontSize = 15;
l = legend('MR MPC','Poly NMPC', 'Nonlinear regulation', 'Feedback linearization');
set(l,'Interpreter','Latex');
l.FontSize = 15;

subplot(3,3,[7,9]);
plot(t*timescale, norm_umpc, 'b','LineWidth', 1.5);
hold on; grid on;
plot(t*timescale, deltaVsr,  'r','LineWidth', 1.5);
plot(t*timescale, norm_ureg,  'color','#D95319', 'LineWidth', 1.5);
plot(t*timescale, norm_ufl,  'g', 'LineWidth', 1.5);
l = xlabel('$Time (h)$');
set(l,'Interpreter','Latex');
l.FontSize = 15;
l = ylabel('Energy expenditure');
set(l,'Interpreter','Latex');
l.FontSize = 15;
% l = legend('MR MPC','Poly NMPC', 'Nonlinear regulation', 'Feedback linearization');
% set(l,'Interpreter','Latex');

