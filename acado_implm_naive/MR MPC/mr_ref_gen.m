function ref = mr_ref_gen(xk,xr,Ts, Np)

ref = zeros(Np,length(xk));
ref1 = xr(1,:)';
ref2 = xr(1,:)';
u = zeros(3,2);
x_t = xk';

delta_b = Ts;
%% First step 

Ad = [ 1, 0, 0, delta_b, 0, 0; ...
      0, 1, 0, 0, delta_b,  0; ...
      0, 0, 1, 0, 0, delta_b; ...
      0, 0, 0, 1, 0,  0; ...
      0, 0, 0, 0, 1, 0; ...
      0, 0, 0, 0, 0, 1];
Am = Ad^2;
  
Bd = [delta_b^2/2 0 0;0 delta_b^2/2 0;0 0 delta_b^2/2; ...
       delta_b 0 0;0 delta_b 0; 0 0 delta_b];

Bm = [Ad*Bd Bd]; 
Cm = [eye(3) zeros(3)];

u1 = pinv(Bm)*(ref1-Am*xk');  
u11 = u1(1); 
u21 = u1(2);
u31 = u1(3);
u12 = u1(4);
u22 = u1(5);
u32 = u1(6);

u_1 = [u11 u12;u21 u22;u31 u32];

%% Second step 
xk1 = Am*xk' + Bm*u1;  
u2 = pinv(Bm)*(ref2-Am*xk1); 

u11 = u2(1); 
u21 = u2(2);
u31 = u2(3);
u12 = u2(4);
u22 = u2(5);
u32 = u2(6);

u_2 = [u11 u12;u21 u22;u31 u32];

%% States planning step
Cd = Cm;
x = x_t;
x1 = Ad*x + Bd*u_1(:,1);
yhat(:,1) = x1;
x2 = Ad*x1 + Bd*u_1(:,2);
yhat(:,2) = x2;
x3 = Ad*x2 + Bd*u_2(:,1);
yhat(:,3) = x3;
x4 = Ad*x3 + Bd*u_2(:,2);
yhat(:,4) = x4;

for i = 0:4:Np
    ref(i+1,:) = yhat(:,1)';
    ref(i+2,:) = yhat(:,2)';
    ref(i+3,:) = yhat(:,3)';
    ref(i+4,:) = yhat(:,4)';
end
end

