% ROBERTSON DAE
% Initialization
mex Jac_create.c -R2018a %compile mex file
pause(0.0001)
clear ; close all; clc
set(0,'DefaultLineLineWidth',2)
warning('off')

%%%%%%%%%%%%%%%%%
%SELECT different tolerances
reltolODE=1e-3;
reltolRPNN=1e-3;
%
abstolODE=reltolODE*1e-3;
abstolRPNN=reltolRPNN*1e-3;

%
t0=0; %initial time
tf=4e6; %final time
u0=[1;0;0]; %initial conditions
%
f=@(t,y) robertson_DAE(t,y); %ODE system
derf=@(t,y) robertson_DAE_jac(t,y); %Jacobian (For RPNN the Jacobian need to be transposed)
%
nEq=max(size(u0)); %number of equations
var_name={'$A$','$B$','$C$'}; %name of variables
%
%Mass Matrix
Mmass=[1,0,0;0,1,0;0,0,0];

%TRUE/REFERENCE solution
iter=1;
opts_true = odeset('Mass',Mmass,'RelTol',1e-14,'AbsTol',1e-16,'Jacobian',derf);
tstart=tic;
sol_true=ode15s(f,[t0,tf],u0,opts_true);
timetrue=toc(tstart);
stepstrue=length(sol_true.x); %number of steps computed
%
tspan=sol_true.x; % the dense output time tspan (used for comparison)
utrue=sol_true.y;
%call RPNN
optsRPNN.Mass=Mmass;
optsRPNN.RelTol=reltolRPNN;
optsRPNN.AbsTol=abstolRPNN;
optsRPNN.Jacobian=derf;
tstart=tic;
[TT,uRPNN,info]=ada_RPNN_DAE(f,tspan,u0,optsRPNN);
timeRPNN=toc(tstart); %execution time
L2errRPNN=norm(uRPNN-utrue,2);
errRPNN=abs(uRPNN-utrue);
stepsRPNN=info.num_steps; %number of steps computed
%
%
%odesuit
opts=odeset('Mass',Mmass,'RelTol',reltolODE,'AbsTol',abstolODE,'Jacobian',derf);
%(ode45 cannot solve DAEs)
%call ode15s
tstart=tic;
sol15s=ode15s(f,[t0,tf],u0,opts);
u15s=deval(sol15s,tspan);
time15s=toc(tstart); %execution times
L2err15s=norm(u15s-utrue,2);
err15s=abs(u15s-utrue);
steps15s=length(sol15s.x); %number of steps computed
%call ode23t
tstart=tic;
sol23t=ode23t(f,[t0,tf],u0,opts);
u23t=deval(sol23t,tspan);
time23t=toc(tstart); %execution times
L2err23t=norm(u23t-utrue,2);
err23t=abs(u23t-utrue);
steps23t=length(sol23t.x); %number of steps computed

            
%FIGURES
for i=1:nEq
figure(i)
plot(tspan,utrue(i,:),'k');
hold on
plot(tspan,u15s(i,:),'--');
plot(tspan,u23t(i,:),'-.');
plot(tspan,uRPNN(i,:),':');
legend('reference','ode15s','ode23t','RPNN')
xlabel('t','interpreter','latex')
ylabel(['$',var_name{i},'$'],'interpreter','latex')
set(gca,'FontSize',16)
%
figure(nEq+i) %errors
semilogy(tspan,err15s(i,:),'--');
hold on
semilogy(tspan,err23t(i,:),'-.');
semilogy(tspan,errRPNN(i,:),':');
legend('ode15s','ode23t','RPNN')
xlabel('t','interpreter','latex')
ylabel(sprintf('$|\\bar{%s}-%s|$',var_name{i},var_name{i}),'interpreter','latex')
set(gca,'FontSize',16)
end

%TABLE
format shorte
method={'ode15s';'ode23t';'RPNN';'reference'};
reltol=[reltolODE;reltolODE;reltolRPNN;1e-14];
L2err={L2err15s;L2err23t;L2errRPNN;0};
exec_time=[time15s;time23t;timeRPNN;timetrue];
steps=uint16([steps15s;steps23t;stepsRPNN;stepstrue]);
T=table(method,reltol,L2err,exec_time,steps)
format short
