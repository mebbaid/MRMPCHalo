clc;
clear all; close all;


%% PARAMETERS
Ts = 0.05;
EXPORT = 1;
mu = 0.012149;
L2 = 1.556;
insertion_error = 0; 
c = zeros(2,1); c(1) = (1-mu)/(L2 + mu)^3 + mu/(L2 - 1 + mu)^3; c(2) = 1/(L2 + mu)^3 - 1/(L2 - 1 + mu)^3;
h = zeros(2,1); h(1) = -1/(2*c(1)*(1-c(1)))*((2-c(1))*4*(L2-2*c(2)*mu*(1-mu))-4*L2);
h(2) = -1/(2*c(1)*(1-c(1)))*(-2*(4*L2-2*c(2)*mu*(1-mu))+4*L2*(1+c(1)));
k = 0.1; 
Omega = 1.8636;
Omega_z = Omega; 
a = zeros(2,1); b = zeros(2,1); a(1) = -1; a(2) = 1/2; b(1) = 2; b(2) = 5/2;
phi = 0; e   = 0; sim_time = 15;

DifferentialState xx xy xz xxd xyd xzd; % Differential States: states
Control ux uy uz; % Control:

%% Differential Equation
d1r = [-mu/(1-mu);0;0];
d2r = [1-mu;0;0];
ed1 = [xx;xy;xz] - d1r;
absed1 = sqrt(ed1(1)^2+ed1(2)^2+ed1(3)^2);
ed2 = [xx;xy;xz] - d2r;
absed2 = sqrt(ed2(1)^2+ed2(2)^2+ed2(3)^2);

f_expl = [ dot(xx) == xxd; ...
           dot(xy) == xyd; ...
           dot(xz) == xzd; ...
           dot(xxd) == 2*xyd+xx-(1-mu)*(xx+mu)/absed1^3 - mu*(xx-1+mu)/absed2^3 + ux; ...
           dot(xyd) == -2*xxd + xy - (1-mu)*xy/absed1^3 - mu*xy/absed2^3 + uy; ...
           dot(xzd) == -(1-mu)*xz/absed1^3 - mu*xz/absed2^3 + uz];% SIMexport
       
acadoSet('problemname', 'sim_crtbp');
numSteps = 3;
sim = acado.SIMexport( Ts );
sim.setModel(f_expl);
sim.set( 'INTEGRATOR_TYPE', 'INT_IRK_GL2' );
sim.set( 'NUM_INTEGRATOR_STEPS', numSteps );
if EXPORT
    sim.exportCode( 'export_SIM' );
    cd export_SIM
    make_acado_integrator('../integrate_crtbp')
    cd ..
end


%% MPCexport
acadoSet('problemname', 'nmpc');
Np = 40;
ocp = acado.OCP( 0.0, Np*Ts, Np );
h = [xx xy xz xxd xyd xzd ux uy uz];
hN = [xx xy xz xxd xyd xzd]; % terminal penalty
W = acado.BMatrix(eye(length(h)));
WN = acado.BMatrix(eye(length(hN))); % terminal penalty weights
ocp.minimizeLSQ( W, h );
ocp.minimizeLSQEndTerm( WN, hN );
ocp.setModel(f_expl);
mpc = acado.OCPexport( ocp );
mpc.set( 'HESSIAN_APPROXIMATION', 'GAUSS_NEWTON' );
mpc.set( 'DISCRETIZATION_TYPE', 'MULTIPLE_SHOOTING' );
mpc.set( 'SPARSE_QP_SOLUTION', 'FULL_CONDENSING_N2');
mpc.set( 'LEVENBERG_MARQUARDT', 1e-5 );
mpc.set( 'INTEGRATOR_TYPE', 'INT_RK4' );
mpc.set( 'NUM_INTEGRATOR_STEPS', numSteps*Np );
mpc.set( 'QP_SOLVER', 'QP_QPOASES' );

if EXPORT
    mpc.exportCode( 'export_MPC' );
    global ACADO_;
    copyfile([ACADO_.pwd '/../../external_packages/qpoases'], 'export_MPC/qpoases')
    cd export_MPC
    make_acado_solver('../acado_MPCstep')
    cd ..
end

%% SIMULATION
% SIM PARAMETERS
X0 = [L2 1 0 0 0 0];
% ref gen
Seq = [L2 + insertion_error;0;0;0;0;0];   
ho1 = [(-k*(1-c(1)+Omega^2)/(2*Omega))*cos(0);
           k*sin(0);
           k*cos(0)];

diffho1 = [(k*(1-c(1)+Omega^2)/2)*sin(Omega*0);
               Omega*k*cos(Omega*0);
               -Omega_z*k*sin(Omega_z*0)];

xr = Seq + [ho1;diffho1];
Xref      = xr';
input.x = [repmat(X0,Np/2,1); repmat(Xref,Np/2+1,1)];
input.od = [];
Uref = zeros(Np,3);
input.u = Uref;
input.y = [repmat(Xref,Np,1) Uref];
% input.y = [xr(:,2:end)' Uref];
input.yN = xr';
input.W = diag([10 10 10 1 1 1 0.1 0.1 0.1]);
input.WN = diag([10 10 10 1 1 1]); 
input.shifting.strategy = 1;

% SIMULATION LOOP
display('------------------------------------------------------------------')
display(' Simulation Loop' )
display('------------------------------------------------------------------')

iter = 0; time = 0;
KKT_MPC = []; ValueFunc_MPC = []; INFO_MPC = [];
controls_MPC = [];
pos_ref = Xref(1:3);
state_sim = X0;
int_step  = 10^-4;

while time(end) < sim_time
    tic
    % Solve NMPC OCP
    input.x0 = state_sim(end,:);
    output = acado_MPCstep(input);
    % Save the MPC step
    INFO_MPC = [INFO_MPC; output.info];
    KKT_MPC = [KKT_MPC; output.info.kktValue];
    ValueFunc_MPC = [ValueFunc_MPC; output.info.objValue];
    QPitr_MPC = [ValueFunc_MPC; output.info.QP_iter];
    controls_MPC = [controls_MPC; output.u(1,:)];
    input.x = output.x;
    input.u = output.u;
    % Simulate system
    sim_input.x = state_sim(end,:).';
    sim_input.u = output.u(1,:).';
    sim_input.od = int_step;
    states = integrate_crtbp(sim_input);
    state_sim = [state_sim; states.value']; 
    iter = iter+1;
    nextTime = iter*Ts;
    disp(['current time: ' num2str(nextTime) ' ' char(9) ' (RTI step -- QP error: ' num2str(output.info.status) ',' ' ' char(2) ' KKT val: ' num2str(output.info.kktValue,'% 1.2e') ',' ' ' char(2) ' CPU time: ' num2str(round(output.info.cpuTime*1e6)) ' µs)'])
    time = [time nextTime];
    % update reference
    ho1 = [(-k*(1-c(1)+Omega^2)/(2*Omega))*cos(nextTime);
           k*sin(nextTime);
           k*cos(nextTime)];

    diffho1 = [(k*(1-c(1)+Omega^2)/2)*sin(Omega*nextTime);
               Omega*k*cos(Omega*nextTime);
               -Omega_z*k*sin(Omega_z*nextTime)];

    xr = Seq + [ho1;diffho1];
    Xref = xr';
    input.y = [repmat(Xref,Np,1) Uref];
    % input.y = [xr(:,2:end)' Uref];
    input.yN = xr';
    pos_ref = [pos_ref ; Xref(1:3)];
end


figure('Name', 'Space craft 3D trajectory');
plot3(state_sim(:,1), state_sim(:,2), state_sim(:,3), 'r','LineWidth' , 2)
hold on; grid on;
plot3(pos_ref(:,1), pos_ref(:,2), pos_ref(:,3), 'b', 'LineWidth' , 2)
scatter3(L2,0,0,'b','diamond');
hold off;

figure;
semilogy(time(1:end-1), ValueFunc_MPC, ':bx');
xlabel('time(s)')
ylabel('Obj_value')
