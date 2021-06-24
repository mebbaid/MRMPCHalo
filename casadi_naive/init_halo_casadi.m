%% Example CASADI MPC for Halo_orbits
% multiple shooting with IP NLP NMPC

clc
clear all


import casadi.*

T  = 0.01; % sampling time
Np = 10 ; % prediction horizon

% parameters
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

Gsc = 1360.8 ; % approx kw
sped = 300000; % approx m/s
zeta = 0.9252 ;   % prep. on the space-craft

u1_min = -inf;  u1_max = inf;
u2_min = -inf;  u2_max = inf;
u3_min = -inf;  u3_max = inf;

% OCP problem parameters
x1 = SX.sym('x1'); x2 = SX.sym('x2'); x3 = SX.sym('x3');
x4 = SX.sym('x4'); x5= SX.sym('x5'); x6 = SX.sym('x6');
states = [x1; x2 ; x3; x4; x5; x6];
n_states = length(states);
u1 = SX.sym('u1'); u2 = SX.sym('u2'); u3 = SX.sym('u3');
controls = [u1; u2 ; u3];
n_controls = length(controls);
z1 = SX.sym('z1'); z2 = SX.sym('z2'); z3 = SX.sym('z3'); z4 = SX.sym('z4'); 
disturbances = [z1;z2;z3;z4]; 
n_dist = length(disturbances);

% dynamics of halo CR3BP prediction model
f1 = [x4;x5;x6];
r1 = sqrt((states(1)+mu)^2+states(2)^2+states(3)^2);
r2 = sqrt((states(1)-1+mu)^2+states(2)^2+states(3)^2);
f2 = [2*states(5)+ states(1) - (1-mu)/r1^3*(states(1)+mu) + mu/r2^3*(states(1)-1+mu); ...
                 -2*states(4) + states(2) - states(2)*(1-mu)/r1^3 - states(2)*mu/r2^3;...
                 -(1-mu)*states(3)/r1^3 - mu*states(3)/r2^3];

rhs  = [f1;f2+controls];

% dynamics of halo ERTBP simulation model

states_1 = states(1:3);
states_2 = states(4:6);
M = [-1 0 0;
    0 -1 0;
    0 0 0];

N = [0 -1 0;
    1 0 0;
    0 0 0];

d1r = [-mu - disturbances(4)/(1-mu);0;0];
d2r = [1-mu + disturbances(4)/mu;0;0];

ed1 = states_1 - d1r;
absed1 = sqrt(ed1(1)^2+ed1(2)^2+ed1(3)^2);
ed2 = states_1 - d2r;
absed2 = sqrt(ed2(1)^2+ed2(2)^2+ed2(3)^2);

f1d = states_2;
f2d = -M*states_1 - 2*N*states_2 - (2*M*states_1 + 2*N*states_2)*disturbances(1) ...
           - M*states_1*disturbances(2) - N*states_1*disturbances(3) - ed1*(1 - mu)/absed1^3 ...
           - ed2*mu/absed2^3;

rhs_d = [f1d;f2d+controls];

f   = Function('f',{states,controls},{rhs});  % nonlinear mapping function object
fd  = Function('fd',{states,controls, disturbances},{rhs_d});  % nonlinear mapping function object

U   = SX.sym('U',n_controls,Np);  % decision variables along horizon
P   = SX.sym('P',n_states+n_states*Np); % parameters including init and ref values of nlp problem
X   = SX.sym('X',n_states,Np+1);     % states along the horizon

% setting OCP and computing obj and equality constraints
obj = 0; % obj function
g   = [];% constraints vector 
Q   = eye(6,6); Q(1,1) = 1000; Q(2,2) = 1000; Q(3,3) = 1000; % weights on states error;
R   = 0.1*eye(3,3); % weights on controls

st  = X(:,1);
g   = [g; st-P(1:n_states)]; % init condition constraints
for k = 1:Np
   st =  X(:,k); con = U(:,k);
   obj = obj+(st-P(n_states*k+1:(k+1)*n_states))'*Q*(st-P(n_states*k+1:(k+1)*n_states)) + con'*R*con;  % quadratic cost along horizon  
   st_next = X(:,k+1);
   f_value = f(st,con);
   st_next_euler = st+(T*f_value);
   g  = [g; st_next - st_next_euler]; % compute path constraints at every point along horizon
end

OPT_variables = [reshape(X,n_states*(Np+1),1) ;reshape(U,n_controls*Np,1)]; % states are also opt variables in multi shooting
nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);  % nlp problem
opts = struct; % solver options
opts.ipopt.max_iter = 1000;
opts.ipopt.print_level = 1;
opts.ipopt.acceptable_tol = 1e-8;

solver = nlpsol('solver', 'ipopt', nlp_prob, opts);

args = struct;  % bounds on the optimization variable
args.lbx(n_states*(Np+1)+1:n_controls:n_states*(Np+1)+n_controls*Np,1) = u1_min; 
args.lbx(n_states*(Np+1)+2:n_controls:n_states*(Np+1)+n_controls*Np,1) = u2_min; % bounds on control along horizon
args.lbx(n_states*(Np+1)+3:n_controls:n_states*(Np+1)+n_controls*Np,1) = u3_min;
args.ubx(n_states*(Np+1)+1:n_controls:n_states*(Np+1)+n_controls*Np,1) = u1_max;  
args.ubx(n_states*(Np+1)+2:n_controls:n_states*(Np+1)+n_controls*Np,1) = u2_max;
args.ubx(n_states*(Np+1)+3:n_controls:n_states*(Np+1)+n_controls*Np,1) = u3_max;

args.lbx(1:n_states:n_states*(Np+1),1) = -inf; % lower bound on x
args.ubx(1:n_states:n_states*(Np+1),1) = inf; % upper bound on x
args.lbx(2:n_states:n_states*(Np+1),1) = -inf; % lower bound on y
args.ubx(2:n_states:n_states*(Np+1),1) = inf; % upper bound on y
args.lbx(3:n_states:n_states*(Np+1),1) = -inf; % lower bound on z
args.ubx(3:n_states:n_states*(Np+1),1) = inf; % upper bound on z
args.lbx(4:n_states:n_states*(Np+1),1) = -inf; % lower bound on xd
args.ubx(4:n_states:n_states*(Np+1),1) = inf; % upper bound on xd
args.lbx(5:n_states:n_states*(Np+1),1) = -inf; % lower bound on yd
args.ubx(5:n_states:n_states*(Np+1),1) = inf; % upper bound on yd
args.lbx(6:n_states:n_states*(Np+1),1) = -inf; % lower bound on zd
args.ubx(6:n_states:n_states*(Np+1),1) = inf; % upper bound on zd


args.lbg(1:n_states*(Np+1)) = 0;   % equality constraints
args.ubg(1:n_states*(Np+1)) = 0;


% Simulation loop
t0 = 0;

Seq = [L2;0;0;0;0;0];   
ho1 = [(-k*(1-c(1)+Omega^2)/(2*Omega))*cos(0);
       k*sin(0);
       k*cos(0)];   
diffho1 = [(k*(1-c(1)+Omega^2)/2)*sin(Omega*0);
           Omega*k*cos(Omega*0);
           -Omega_z*k*sin(Omega_z*0)];      
% x0 = Seq + [ho1;diffho1];
x0 = [L2;0;0;0;0;0]; % ref posture

xx(:,1) = x0;  % history of states
t(1) = t0; % time 
u0 = zeros(Np,n_controls); % two controls along horizon
X0 = repmat(x0,1,Np+1); % init of state decision variables
d0 = [0;0;0;0];


simTime = 15;
mpcitr = 0; % start mpc
xx1 = [];
u_cl = [];
n_ref = [];
dstr = [];

j = 0.1;
% main simulation loop
while(mpcitr < simTime/T)
   args.p(1:n_states) = x0; %values of parameters of opt problem at begining
   current_time = mpcitr*T;
   for k = 1:Np
      t_predict = current_time + (k-1)*T; % predicted time instant
      ho1 = [(-j*(1-c(1)+Omega^2)/(2*Omega))*cos(t_predict);
               j*sin(t_predict);
               j*cos(t_predict)];
      diffho1 = [(j*(1-c(1)+Omega^2)/2)*sin(Omega*t_predict);
               Omega*j*cos(Omega*t_predict);
               -Omega_z*j*sin(Omega_z*t_predict)];
      ref = Seq + [ho1;diffho1];      
      xref = ref(1); yref = ref(2); zref = ref(3);
      xdref = ref(4); ydref = ref(5); zdref = ref(6);
      args.p(n_states*k+1:n_states*(k+1)) = [xref, yref, zref, xdref, ydref, zdref];
   end   
   
   n_ref = [n_ref args.p(n_states*k+1:n_states*(k+1))'];
   args.x0 = [reshape(X0',n_states*(Np+1),1);reshape(u0',n_controls*Np,1)]; % init value of opt variable
   sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx, 'lbg', args.lbg, ...
             'ubg', args.ubg, 'p', args.p);
   u   = reshape(full(sol.x(n_states*(Np+1)+1:end))',n_controls,Np)'; % opt controls in 2*Np 
   xx1(:,1:n_states,mpcitr+1) = reshape(full(sol.x(1:n_states*(Np+1)))',n_states,Np+1)';
   u_cl = [u_cl;u(1,:)];
   t(mpcitr+1) = t0;
   d0 = dist(current_time, a,e, phi, b, mu) ;
   dstr = [dstr d0];
   [t0, x0, u0] = shift(T,t0, x0, u,d0, fd);  % apply control and update opt problem init values x0,u0,t0
   xx(:,mpcitr+2) = x0;
   X0 = reshape(full(sol.x(1:n_states*(Np+1)))',n_states,Np+1)';
   X0 = [X0(2:end,:);X0(end,:)];
   mpcitr = mpcitr + 1;   
end

x = xx(1,:)'; refx = n_ref(1,:)';
y = xx(2,:)'; refy = n_ref(2,:)';
z = xx(3,:)'; refz = n_ref(3,:)';



figure('Name','Trajectories');
subplot(2,1,1);
plot3(refx,refy, refz, 'k', 'LineWidth', 2);
hold on; grid on;
plot3(x,y, z, 'r', 'LineWidth', 2);
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


function [t0, x0,u0] = shift(T,t0, x0, u,d0,fd)
st = x0;
cn = u(1,:)';
dst = d0;
f_value = fd(st,cn,dst);
st = st+ (T*f_value);
x0 = full(st);
t0 = t0 + T;
u0 = [u(2:size(u,1),:); u(size(u,1),:)];
end

function d0 = dist(t, a,e, phi, b, mu)
   a1 = a(1)*e*cos(t+phi);
   alpha = a1;
   b1 = b(1)*e*cos(t+phi);
   bd = -b(1)*e*sin(t+phi);
   beta = [b1;0;bd];
   d0 = [beta;mu*(1-mu)*alpha];
end