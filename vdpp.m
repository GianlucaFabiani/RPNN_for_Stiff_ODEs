function dydt = vdpp(~,y,mu)
%van der Pol ODEs
dydt(1,:)=y(2,:);
dydt(2,:)=-y(1,:)+mu*(1-y(1,:).^2).*y(2,:);
end

