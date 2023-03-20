function J = LotkaVolterra_jac(~,y,U,R)
[~,np]=size(y);
if np==1
    J=[R/U*(2-0.04*U^2*y(2)),R/U*(-0.04*U^2*y(1));...
        R/U*(0.02*U^2*y(2)),R/U*(0.02*U^2*y(1)-1.06*U)];
else
    uu=ones(1,np);
    J=[R/U*(2*uu-0.04*U^2*y(2,:)), R/U*(0.02*U^2*y(2,:));...
        R/U*(-0.04*U^2*y(1,:)),R/U*(0.02*U^2*y(1,:)-1.06*U*uu)];
end
end