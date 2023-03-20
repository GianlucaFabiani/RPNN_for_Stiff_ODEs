function f = robertson_DAE(~,y)
k=[0.04,1e4,3e7]; %parameters of ROBERTSON system
f(1,:)=-k(1)*y(1,:) + k(2)*y(2,:).*y(3,:);
f(2,:)=k(1)*y(1,:) - k(2)*y(2,:).*y(3,:) - k(3)*y(2,:).^2;
%y(3,:)=                                   k(3)*psi(2,:).^2;
f(3,:)=y(1,:) + y(2,:) + y(3,:) - 1; %algebraic for conservation law 1=y1+y2+y3;
end