function J = robertson_DAE_jac(~,y)
k=[0.04,1e4,3e7]; %parameters of ROBERTSON system
[~,np]=size(y);
if np==1
    J=[-k(1),k(2)*y(3),k(2)*y(2);
        k(1),-k(2)*y(3)-2*k(3)*y(2),-k(2)*y(2);...
        1, 1, 1];
else
    uu=ones(1,np);
    J=[-k(1)*uu, k(1)*uu, uu;...
     k(2)*y(3,:), - k(2)*y(3,:) - 2*k(3)*y(2,:), uu;...
      k(2)*y(2,:), -k(2)*y(2,:), uu];
end
end