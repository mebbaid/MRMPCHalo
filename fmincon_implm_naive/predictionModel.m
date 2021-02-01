function xpred = predictionModel(u,x, np,delta)
%SATELLITE Summary of this function goes here
%   Detailed explanation goes here
    ct = u(:,1);
    mu = 0.012149;
    z = zeros(4,1);
        M = [-1 0 0;
        0 -1 0;
        0 0 0];

    N = [0 -1 0;
        1 0 0;
        0 0 0];
    
    x1 = x(1:3);
    x2 = x(4:6);
    xpred = sym('xpred',[6,np],'real');
    xpred(:,1) = x;    
    for j=2:np  
        d1r = [-mu - z(4)/(1-mu);0;0];
        d2r = [1-mu + z(4)/mu;0;0];

        ed1 = x1 - d1r;
        absed1 = sqrt(ed1(1)^2+ed1(2)^2+ed1(3)^2);
        ed2 = x1 - d2r;
        absed2 = sqrt(ed2(1)^2+ed2(2)^2+ed2(3)^2);

        f1 = x2;
        f2 = -M*x1 - 2*N*x2 - (2*M*x1 + 2*N*x2)*z(1) ...
                   - M*x1*z(2) - N*x1*z(3) - ed1*(1 - mu)/absed1^3 ...
                   - ed2*mu/absed2^3;
        F = [f1; f2];
        G = [zeros(3); eye(3)];
        xpred(:,j) = xpred(:,j-1) + delta*( F + G*ct);
        x1 = xpred(1:3,j);
        x2 = xpred(4:6,j);
        ct = u(:,j);
    end 
end


