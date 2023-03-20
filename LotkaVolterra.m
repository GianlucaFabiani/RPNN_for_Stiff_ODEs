function dydt = LotkaVolterra(~,y,U,R)
%U=200;R=20;
dydt(1,:)=R/U*(2*U*y(1,:)-0.04*U^2*y(1,:).*y(2,:));
dydt(2,:)=R/U*(0.02*U^2*y(1,:).*y(2,:)-1.06*U*y(2,:));
end

