function out = ex11ode(t,u,flag)
% Copyright 1999, The MathWorks, Inc.
if nargin < 3 || isempty(flag)
  s = sin(t+pi/4);
  c = cos(t+pi/4);
  out = zeros(5,length(t));
  out(1,:) =   u(2,:);
  out(2,:) = - 10*u(2,:) + u(5,:).*s;
  out(3,:) =   u(4,:);
  out(4,:) = - 10*u(4,:) - u(5,:).*c + 1;
  g     =  c.*u(3,:) - s.*u(1,:);
  gp    =  c.*(u(4,:) - u(1,:)) + s.*( - u(2,:) - u(3,:));
  gpp   =  c.*(out(4,:) - out(1,:) - u(2,:) - u(3,:)) ...
         + s.*(- out(2,:) - out(3,:) - u(4,:) + u(1,:));
  out(5,:) = gpp + 20*gp + 100*g;
else
  switch(flag)           
  case 'mass'
    out = diag([1 1 1 1 0]);
  otherwise
    error(['Unknown flag ''' flag '''.']);
  end
end     