function J=ex11ode_jac(t,y)
s = sin(t+pi/4);
c = cos(t+pi/4);
np=size(y,2);
if np==1
    J=zeros(5,5);
    J(1,2)=1;
    J(2,2)=-10; J(2,5)=s;
    J(3,4)=1;
    J(4,4)=-10; J(4,5)=-c;
    J(5,1) = -99*s-20*c;
    J(5,2) = -10*s-2*c;
    J(5,3) = 99*c-20*s;
    J(5,4) = 10*c-2*s;
    J(5,5) = -c.^2-s.^2; 
else
    ll=length(t);
    J=zeros(5,5*ll);
    ind=1:ll;
    J(2,ind) = 1;

    ind=ll+ind;
    J(2,ind) = -10;
    J(5,ind) = s;

    ind=ll+ind;
    J(4,ind) = 1;

    ind=ll+ind;
    J(4,ind) = -10;
    J(5,ind) = -c;
    ind=ll+ind;
    J(1,ind) = -99*s-20*c;
    J(2,ind) = -10*s-2*c;
    J(3,ind) = 99*c-20*s;
    J(4,ind) = 10*c-2*s;
    J(5,ind) = -c.^2-s.^2; 
end
end