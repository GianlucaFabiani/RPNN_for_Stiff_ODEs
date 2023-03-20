function [tt_tot,uu_tot,info]=ada_RPNN_DAE(f,tspan,y0,opts)
% Physics-Informed Random Projection Neural Network for solving system of
% stiff ODEs and DAEs
%
% G. Fabiani, March 15, 2023
%
% [TOUT,YOUT] = ada_RPNN_DAE(ODEFUN,TSPAN,Y0)
%
% integrates the system of differential equations
% y' = f(t,y), from time T0 to TFINAL and y(T0)=Y0.
%
% INPUTS:
% ODEFUN: function handle
% For a row-vector of times T and a matrix Y (rows: times,colums:
% equations) ODEFUN(T,Y) must return a matrix, with columns
% corresponding to f(t,y).
% - TSPAN = [T0 TFINAL]
% To obtain solutions at specific times T0,T1,...,TFINAL (all increasing), use TSPAN =
% [T0 T1 ... TFINAL].
% - Y0: initial conditions, should include also the starting value for DAEs
%
% OUTPUTS: Each column in the solution array YOUT corresponds
% to a time returned in the row vector TOUT.
%
% [TOUT,YOUT] = ode15s(ODEFUN,TSPAN,Y0,OPTIONS) solves as above with default
% integration properties replaced by values in OPTIONS (struct with fields)
% Commonly used options are relative error tolerance 'RelTol'
% (1e-3 by default) and absolute error tolerances 'AbsTol' (1e-6 by default);
%
% As output one can also get a third argument for info about the computations
% [TOUT,YOUT,INFO] = ode15s(ODEFUN,TSPAN,Y0,OPTIONS)
% INFO.accept: number of accepted adaptive steps
% INFO.reject: number of rejected adaptive steps
% INFO.exec_time: execution time
% INFO.num_steps: number of total step needed
% INFO.num_fun_eval: number of function evaluations
% INFO.num_Jac_solv: number of jacobian inversions(pinv/QR decompositions)

if length(tspan)==2
    flagtspan=0; %case tspan=[t0,tf]; The output is written in adaptive found time points.
    if length(tspan)==1
        tspan=[0,tspan];
    end
else
    flagtspan=1; %case tspan=[t0,t1,...,tf], where the output is sought.
end
tf=tspan(end);
t0=tspan(1);
if nargin==3
 opts=[]; %default options
end
nEq=max(size(y0)); %number of equations

if ~isfield(opts,'AbsTol') %Absolute Tolerance
    abstol=1e-6;
else
    if isempty(opts.AbsTol)
        abstol=1e-6;
    else
        abstol=opts.AbsTol;
    end
end
if ~isfield(opts,'RelTol') %Relative Tolerance
    reltol=1e-3;
else
    if isempty(opts.RelTol)
        reltol=1e-3;
    else
        reltol=opts.RelTol;
    end
end
if ~isfield(opts,'Jacobian')
    dfdu='FD'; %the Jacobian is approximated with finite differences
else
    if isempty(opts.Jacobian)
        dfdu='FD';
    else
        dfdu=opts.Jacobian; %analytical jacobian is provided
% for the current version it needs function handle that gives in output
% the transposed Jacobian that can be evaluated in more than one point in time
    end
end
if ~isfield(opts,'N') %number of Neurons (the code is optimized only for N=20)
    N=20;
else
    if isempty(opts.N)
        N=20;
    else
        N=opts.N;
    end
end
if ~isfield(opts,'np') %number of collocation points (the code is optimized only for np=20)
    nnt=20;
else
    if isempty(opts.np)
        nnt=20;
    else
        nnt=opts.np;
    end
end
if ~isfield(opts,'Method') %selection of the method for the solution of the
    %linear system in each iteration of the Newton scheme
    if nnt*N*nEq^2<20000
        method=1; %pinv or decomposition
    else
        method=2; %sparse qr
    end
else
    if isempty(opts.Method)
        if nnt*N*nEq^2<20000
            method=1; %pinv or decomposition
        else
            method=2; %sparse qr
        end
    else
        method=opts.Method;
    end
end
if ~isfield(opts,'Mass') %Mass Matrix for DAEs
    opts.Mass=[];
    if method==2 %sparse
        Mmass=speye(nEq);
        [Mmassi,Mmassj,Mmassval]=find(Mmass);
    else %full
        Mmass=eye(nEq);
    end
else
    if isempty(opts.Mass)
        if method==2
            Mmass=speye(nEq);
            [Mmassi,Mmassj,Mmassval]=find(Mmass);
        else
            Mmass=eye(nEq);
        end
    else
        Mmass=opts.Mass;
        if method==2 && ~issparse(Mmass)
            Mmass=sparse(Mmass);
        end
        if method==2
            [Mmassi,Mmassj,Mmassval]=find(Mmass);
        end
    end
end
%activation function of the network (optimized only for default Gaussian rbf)
tr_fun=@rbf;
tr_fun_t=@rbf_t;
valw=1;

if isempty(opts.Mass)
    % function for initializing first step usally used for Dormand-Prince
    % RK5(4)
    DT=initialize_first_step(f,y0,t0,abstol,reltol)*10;
    dy1dt=f(t0,y0); %initial slope
else
    % modification for DAEs (not optimal)
    [DT,dy1dt]=initialize_first_stepDAE(Mmass,f,y0,t0,abstol,reltol);
    DT=DT*10;
end

%adaptation parameters
facmin=0.1;
fac=0.8;
facmax=5;
reject=0;
accept=0;
der_prec=dy1dt;

cen_0=linspace(0,1,N)'; %equispace centers (of the N Gaussians)
tk1=t0;
tt_tot=[];
uu_tot=[];
IC=y0; %initial condition (that will change in each subinterval)
ks=0; %number of adaptation
maxiter=5; %total maxiter
if method==1
    %number of true Jacobian computed
    maxiter2=max(min(floor(3-1/3-log10(reltol)/3),maxiter),1);
elseif method==2
    maxiter2=maxiter; %for the suiteSparse is not optimal quasi-newton
%as it seems to be implemented
end

%these go in output as "info"
num_steps=0; %number of adaptive steps
num_fun_eval=0; %number of function evaluations
num_Jac_solv=0; %number of Jacobian pinv/QR decomposition
time_RPN=0; %time for only solving (not evaluating the solution in new points)
while tf-tk1>100*eps(tk1)
    tstart=tic;
    % Stretch the step if within 10% of tfinal-t.
    if 1.1*DT>(tf-tk1)
        DT=tf-tk1;
    end
    ks=ks+1;
    err=100;
    while err>1
        tk2=tk1+DT;
        tt=linspace(tk1,tk2,nnt);
        t_tk=tt-tk1;
        %Parsimonius RBF
        centers=cen_0*DT+tk1;      % centers of RBFs
        amax=25/(9*DT^2); %best bound for the uniform distribution
        s=amax*rand(N,1);
        %output hidden layer
        PHI=tr_fun(tt,centers,s);
        PHI_t=tr_fun_t(tt,centers,s);
        PHIder=PHI+t_tk.*PHI_t;
        PHIt_tk=-t_tk.*PHI;
        
        %initial guess
        if any(der_prec)
            PHIend_temp=PHI(:,end);
            w=der_prec.*(PHIend_temp'/(sum(PHIend_temp.^2))); %continuation
            % in the parameter space
            % looking for u(Tfinal)=u(T0)+dt*u_t(T0)
            % and comparing with
            % u(Tfinal)=u(T0)+dt*w*PHI(Tfinal)
        else
            w=zeros(nEq,N);
        end
        out=w*PHI;
        out_t=w*PHI_t;
        psi=IC+t_tk.*out; %trial solution
        psi_t=out+t_tk.*out_t; %time derivative
        eqn=f(tt,psi); %rhs of equations
        num_fun_eval=num_fun_eval+length(tt);
        if isempty(opts.Mass)
            G=eqn-psi_t;
        else
            G=(eqn-Mmass*psi_t); %G=Mmass*u_t-eqn, -G.
        end
        Gre=reshape(G',1,[]);
        %err=norm(Gre)
        %Newton method
        errNew=100;
        if isempty(opts.Mass)
            err=norm( Gre./ (abstol+reltol*reshape(abs(psi_t)',1,[])) );
        else
            err=norm( Gre./ (abstol+reltol*reshape(abs(Mmass*psi_t)',1,[])) );
        end
        err_prec=1/abstol*1e5;
        iter=0;
        %NO gmres
        while iter<maxiter && err_prec>err && (errNew>1e-14 && err>1)
            err_prec=err;
            iter=iter+1;
            if iter<=maxiter2
                if ischar(dfdu)
                    %ep=1e-9;
                    psi_p=psi;
                    %psi_m=psi;
                    DD=zeros(nEq,nEq*nnt);
                    for i=1:nEq
                        ep=sqrt(eps(psi(i,:)));
                        psi_p(i,:)=psi(i,:)+ep;
                        %psi_m(i,:)=psi(i,:)-ep;
                        ffp=f(tt,psi_p)./(ep);
                        num_fun_eval=num_fun_eval+length(tt);
                        ffa=eqn./ep;
                        %ffm=reshape((f(tt,psi_m)./(2*ep))',1,[]);
                        %DD(i,:)=(ffp-ffm);
                        DD(i,:)=reshape((ffp-ffa)',1,[]);
                        psi_p(i,:)=psi(i,:);
                        %psi_m(i,:)=psi(i,:);
                    end
                else
                    DD=dfdu(tt,psi);
                end
                if method==2
                    if ~issparse(DD)
                        DD=sparse(DD);
                    end
                    [DDi, DDj, DDval]=find(DD);
                    [J1i,J1j,J1val,J2i,J2j,J2val] = Jac_create_sparse(Mmassi,Mmassj,Mmassval, PHIt_tk, DDi, DDj, DDval, PHIder);
                    J1=sparse(J1i,J1j,J1val,nEq*N,nEq*nnt);
                    J2=sparse(J2i,J2j,J2val,nEq*N,nEq*nnt);
                    J=J1+J2;
                else
                    method=1;
                    if issparse(DD)
                        DD=full(DD);
                    end
                    J=Jac_create(Mmass,PHIt_tk,DD,PHIder);
                end
                t1=0;
                if method==1
                    if nEq<8
                        invJ=pinv(J);
                        num_Jac_solv=num_Jac_solv+1;
                    else
                        dA=decomposition(J,'cod');
                        num_Jac_solv=num_Jac_solv+1;
                    end
                end
            end
            if method==2
                dwre=spqr_solve(J',Gre', struct('solution','min2norm'));
                num_Jac_solv=num_Jac_solv+1;
            elseif method==1
                if nEq<8
                    dwre=Gre*invJ; %implicit -G (G*=eqn-u_t=-G)
                else
                    dwre=Gre/dA;
                end
            elseif method==0
                [dwre,~]=gmres(J',Gre',[],1e-10,4*ceil(sqrt(nEq*N)));
                num_Jac_solv=num_Jac_solv+1;
            end
            dw=reshape(dwre,[],nEq)';
            w=w+dw;
            errNew=norm(t_tk.*(dw*PHI));
            out=w*PHI;
            out_t=w*PHI_t;
            psi=IC+t_tk.*out;
            psi_t=out+t_tk.*out_t;
            eqn=f(tt,psi);
            num_fun_eval=num_fun_eval+length(tt);
            if isempty(opts.Mass)
                G=(eqn-psi_t);
            else
                G=(eqn-Mmass*psi_t); %f(t,u) dfdw=df/du*du/dw
            end
            Gre=reshape(G',1,[]);
            if isempty(opts.Mass)
                err=norm( Gre./ (abstol+reltol*reshape(abs(psi_t)',1,[])) );
            else
                err=norm( Gre./ (abstol+reltol*reshape(abs(Mmass*psi_t)',1,[])) );
            end
            %fprintf('iter=%d, errG=%2.4e, t=%2.6f, N=%d, dt=%2.4e\n',iter,err,tk1,N,DT);
        end
        if err<1
            reject=reject+1;
        end
        par=(1/(err))^(1/(iter+1));
        if par>1 && iter<maxiter2
            par=par*(maxiter2/iter);
        end
        step_inc=fac*min(facmax,max(facmin,par));
        DT=DT*step_inc;
    end
    num_steps=num_steps+1;
    %output end
    PHIend=tr_fun(tk2,centers,s);
    PHIend_t=tr_fun_t(tk2,centers,s);
    %update IC
    ICprev=IC;
    IC=IC+(tk2-tk1).*(w*PHIend);
    %der_prec=f(tend,IC);
    der_prec=w*PHIend+(tk2-tk1).*(w*PHIend_t);
    der_prec(abs(der_prec)<abstol^2+1e-10)=0;
    accept=accept+1;
    time_el=toc(tstart);
    time_RPN=time_RPN+time_el;
    
    if flagtspan==1
    t_tot=tspan(tspan>=tk1);% && tspan<tk2));
    t_tot=t_tot(t_tot<tk2);
    if ~isempty(t_tot)
        if tf-tk2<= 100*eps(tk2) && t_tot(end)<tf
            t_tot=[t_tot,tf];
        end
    elseif tf-tk2<= 100*eps(tk2)
        t_tot=tf;
    end
    tt_tot=[tt_tot,t_tot]; %%%
    %total coordinates
    if ~isempty(t_tot)
        %output tot
        PHItot=tr_fun(t_tot,centers,s);
        %update u_tot
        uu=ICprev+(t_tot-tk1).*(w*PHItot);
        uu_tot=[uu_tot,uu];
    end
    elseif flagtspan==0
        if ~isempty(tt_tot)
            tt=tt(2:end);
            psi=psi(:,2:end);
        end
        tt_tot=[tt_tot,tt];
        uu_tot=[uu_tot,psi];
    end
    tk1=tk2;
end
info=[];
info.accept=accept;
info.reject=reject;
info.exec_time=time_RPN;
info.num_steps=num_steps;
info.num_fun_eval=num_fun_eval;
info.num_Jac_solv=num_Jac_solv;
end

function g = rbf(t,c,s)
g = exp(-s.*((t-c).^2));
end

function g = rbf_t(t,c,s)
g = rbf(t,c,s).*(-2*s.*(t-c));
end

function [h,dy1dt]=initialize_first_stepDAE(Mmass,fname,y0,t0,AbsTol,RelTol)
nEq=max(size(y0));
% Inizialization of the step h
sc=AbsTol+RelTol*y0;
d0=1/nEq*norm(y0./sc);
ff=fname(t0,y0);
if issparse(Mmass)
    Mmass=full(Mmass);
end
PP=pinv(Mmass);
hh=1e-8;
ep=1e-8;
t1=t0+hh;
y1=y0+hh*PP*ff;
ff1=fname(t1,y1);
y1=y0+hh*PP*(ff+ff1)/2;
ff1=fname(t1,y1);
for kk=1:5
F=Mmass*(y1-y0)-hh*(ff+ff1)/2;
if kk==1
y1p=y1;
y1m=y1;
J=zeros(nEq,nEq);
for i=1:nEq
    y1p(i)=y1(i)+ep;
    y1m(i)=y1(i)-ep;
    ff1p=fname(t1,y1p);
    ff1m=fname(t1,y1m);
    J(:,i)=(ff1p-ff1m)/(4*ep);
    y1p(i)=y1(i);
    y1m(i)=y1(i);
end
PJ=pinv(Mmass-hh*J);
end
dy1=-PJ*F;
y1=y1+dy1;
ff1=fname(t1,y1);
end
dy1dt=(y1-y0)/hh;
d1=1/nEq*norm(dy1dt./sc);
if (min(d0,d1)<10^(-5)) %check if d0 or d1 are too small
    h0=10^(-6);
else
    h0=0.01*d0/d1;  %first hypotesis on step h
end

%Forward Euler method with step h0
y1=y0+h0*PP*ff;
t1=t0+h0;
ff1=fname(t1,y1);
y1=y0+h0*PP*(ff+ff1)/2;
for kk=1:3
F=Mmass*(y1-y0)-h0*(ff+ff1)/2;
if kk==1
y1p=y1;
y1m=y1;
for i=1:nEq
    y1p(i)=y1(i)+ep;
    y1m(i)=y1(i)-ep;
    ff1p=fname(t1,y1p);
    ff1m=fname(t1,y1m);
    J(:,i)=(ff1p-ff1m)/(4*ep);
    y1p(i)=y1(i);
    y1m(i)=y1(i);
end
PJ=pinv(Mmass-h0*J);
end
dy1=-PJ*F;
y1=y1+dy1;
ff1=fname(t1,y1);
end
%%%%%%%
hh=1e-8;
ep=1e-8;
t2=t1+hh;
y2=y1+hh*PP*ff1;
ff2=fname(t2,y2);
y2=y1+hh*PP*(ff1+ff2)/2;
ff2=fname(t2,y2);
for kk=1:3
F=Mmass*(y2-y1)-hh*(ff1+ff2)/2;
if kk==1
y2p=y2;
y2m=y2;
for i=1:nEq
    y2p(i)=y2(i)+ep;
    y2m(i)=y2(i)-ep;
    ff2p=fname(t2,y2p);
    ff2m=fname(t2,y2m);
    J(:,i)=(ff2p-ff2m)/(4*ep);
    y2p(i)=y2(i);
    y2m(i)=y2(i);
end
PJ=pinv(Mmass-hh*J);
end
dy2=-PJ*F;
y2=y2+dy2;
ff2=fname(t2,y2);
end
dy2dt=(y2-y1)/hh;

d2=1/nEq*norm(dy2dt./sc);
d2=d2/h0;
if max(d1,d2)<10^(-15) %check d1 and d2
    h1=max(10^(-6),h0*10^(-3));
else
    h1=0.01/max(d1,d2);
    h1=h1^(1/6);  %second hypothesis on step h
end
h=min(100*h0,h1); %computed initial step h
%h=5*h;
end

function [h]=initialize_first_step(fname,y0,t0,AbsTol,RelTol)
m=max(size(y0));
% Inizialization of the step h
sc=AbsTol+RelTol*y0;
d0=1/m*norm(y0./sc);
ff=fname(t0,y0);
d1=1/m*norm(ff./sc);
if (min(d0,d1)<10^(-5)) %check if d0 or d1 are too small
    h0=10^(-6);
else
    h0=0.01*d0/d1;  %first hypotesis on step h
end

%Forward Euler method with step h0
y1=y0+h0*fname(t0,y0);
f1=fname(t0+h0,y1);
d2=1/m*norm(f1./sc);
d2=d2/h0;
if max(d1,d2)<10^(-15) %check d1 and d2
    h1=max(10^(-6),h0*10^(-3));
else
    h1=0.01/max(d1,d2);
    h1=h1^(1/6);  %second hypothesis on step h
end
h=min(100*h0,h1); %computed initial step h
%h=5*h;
end
