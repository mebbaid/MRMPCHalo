%% main loop preparing ode and OCP.
clc; clear all;

% system parameters
L2      = 1.1556;
mu      = 0.012149;
c       = zeros(2,1);
c(1)    = (1-mu)/(L2 + mu)^3 + mu/(L2 - 1 + mu)^3;
c(2)    = 1/(L2 + mu)^3 - 1/(L2 - 1 + mu)^3;
h       = zeros(2,1);
h(1)    = -1/(2*c(1)*(1-c(1)))*((2-c(1))*4*(L2-2*c(2)*mu*(1-mu))-4*L2);
h(2)    = -1/(2*c(1)*(1-c(1)))*(-2*(4*L2-2*c(2)*mu*(1-mu))+4*L2*(1+c(1)));
k       = 0.1; 
Omega   = 1.8636;
Omega_z = Omega; 
a       = zeros(2,1);
b       = zeros(2,1);
a(1)    = -1;
a(2)    = 1/2;
b(1)    = 2;
b(2)    = 5/2;
phi     = 0;

insertion_error = 0;  % orbit insertion error (set to -0.03 to start near orbit)
distanceScale = 384400; % distance between two primaries
% errorScale = distanceScale/1000;  % error in km

% simulation parameters
disturbance = 0;
saturation  = 0;
sat_constraint = 0;
if saturation == 1
    satValue = 15; % adjust this value to saturate inputs
else
    satValue = inf;
end

t0       = 0;
delta    = 0.05; % sampling period
np       = 4;   % prediction horizon
simTime  = 10;
planttol = 10^-5;
predtol  = delta;
nx       = 6;
ny       = nx;
nu       = 3;
nd       = 4;
Q        = diag([10 10 10 1 1 1]);
R        = 0.1*eye(nu);


x0      = [L2;0;0;0;0;0];
u0      = zeros(nu,np);
cost    = 0;
% xpred   = zeros(nx,np);
xr      = zeros(nx,1);
% xpred(:,1) = x0;
options    = optimoptions('fmincon','Display','iter','Algorithm','sqp');
odeoptions = odeset('RelTol',10^-5, 'AbsTol',1e-7);

% nonlinear constraints
c_ineq = [];
c_eq   = [];

% linear constraints
A = []; Aeq = [];
B = []; Beq = [];

% bounds
if sat_constraint == 1
    lb = -satValue*ones(nu,np);
    ub = satValue*ones(nu,np);
else
    lb = [];
    ub = [];
end

% ref preparation
Seq     = [L2 + insertion_error;0;0;0;0;0];   
e       = 0.0549;
ho1     = zeros(3,np);
diffho1 = zeros(3,np);

% plant and disturbances preparation
u = sym('u',[nu np], 'real');
x = zeros(nx,1);
z = zeros(nd,1);

% save data for plotting later
% save data for plotting later
xtraj = x0;
ztraj = [];
utraj = [];
erms  = [];
deltaVpoly = [];
ref   = [];
PolyNmpciter = [];


% main loop
for i=0:delta:simTime
   xpred =predictionModel(u,x0,np,delta);
   xr = ref_gen(xr,np,k,c,Omega,Omega_z,delta,Seq,i);
   cost = costFunc(u,xpred,[xr xr xr xr],Q,R,np);
   fun = matlabFunction(cost,'Vars',{u});
   % prepare and solve OCP
   [ct,fval,exitflag,fminconoutput] = fmincon(fun,u0,A,B,Aeq,Beq,lb,ub,[],options);
   % apply control, update states
   if disturbance ==1
        z = perturbation(phi,i, e,a,b,mu);
   end
   if saturation ==1
        for j = 1:nu
            ct(j,1) = min(satValue, max(ct(j,1), -satValue));   % implements saturation with script
        end
   end
   [s,x] = ode45(@(t,x) plant(i,x,z,ct(:,1)),[i, i+delta],x0,odeoptions);
   % update and store data
   x = x(length(x),:)';
   er = x(1:3)-xr(1:3);
   xpred(:,1) = x;
   xtraj = [xtraj x]; ztraj = [ztraj z];
   deltaVpoly = [deltaVpoly sqrt(ct(:,1)'*ct(:,1))];
   utraj = [utraj ct(:,1)];
   erms = [erms sqrt(er'*er)];
   ref   = [ref xr];
   xr    = zeros(nx,1);
   x0         = x;
   PolyNmpciter = [PolyNmpciter fminconoutput.iterations];
   u0         = ct;  %not always wise to initialize fmincon with prev control, might get stuck on unfeasible solution
end


%% plotting
t =0:delta:simTime;

figure('Name','MPC using polyn, pred. model');

subplot(2,2,1);
l = title('Three dimensional plot');
set(l,'Interpreter','Latex');
plot3(xtraj(1,:),xtraj(2,:), xtraj(3,:), 'r', 'LineWidth', 1.5);
hold on; grid on;
plot3(ref(1,:),ref(2,:), ref(3,:), 'k', 'LineWidth', 1.5);
hold on
scatter3(L2,0,0,'b','diamond');
l = xlabel('$x_1$');
set(l,'Interpreter','Latex');
l = ylabel('$x_2$');
set(l,'Interpreter','Latex');
l = zlabel('$x_3$');
set(l,'Interpreter','Latex');
l = legend('$x(t), y(t), z(t)$- NMPC trajectory','$x(t), y(t), z(t)$- Nominal trajectory', 'L2 point' );


subplot(2,2,2)
l = title('Controls');
set(l,'Interpreter','Latex');
plot(t, utraj(1,:), 'k', 'LineWidth', 1.5);
hold on; grid on;
plot(t, utraj(2,:), 'r', 'LineWidth', 1.5);
plot(t, utraj(3,:), 'b', 'LineWidth', 1.5);
l = legend('$u_1(t)$', '$u_2(t)$', '$u_3(t)$');
set(l,'Interpreter','Latex');
l = xlabel('Time (h)'); 
l.FontSize = 18;


subplot(2,2,3);
l = title('Norm of the error');
set(l,'Interpreter','Latex');
plot(t, erms, 'r', 'LineWidth', 1.5);
hold on; grid on;
l = legend('NMPC $\|e(t)\|$');
set(l,'Interpreter','Latex');
l = xlabel('Time (h)'); 
l.FontSize = 18;
hold off;


subplot(2,2,4);
l = title('Control effort magnitutde');
set(l,'Interpreter','Latex');
plot(t, deltaVpoly, 'k', 'LineWidth', 1.5);
hold on; grid on;
l = legend('NMPC fmincon $\|u\|(t)$');
set(l,'Interpreter','Latex');
l = xlabel('Time (h)'); 
l.FontSize = 18;
hold off;