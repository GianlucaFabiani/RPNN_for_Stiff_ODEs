function J = vdp_jac(~,y,mu)
[~,np]=size(y);
if np==1
    J=[0 , 1;...
    -1-2*mu*y(1)*y(2), mu*(1-y(1)^2)];
else % transposed jacobian evaluated in np points
    zz=zeros(1,np);
    one=ones(1,np);
    J=[zz, -one-2*mu*y(1,:).*y(2,:);...
    one,mu*(1-y(1,:).^2)];
end
end