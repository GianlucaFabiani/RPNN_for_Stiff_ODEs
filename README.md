# RanDiffNet MATLAB TOOLBOX
# Random Projection PINNs for the Numerical Solution of Stiff ODEs and DAEs

*If you use or modify for research purposes this software, please cite our paper as below:*

&#x1F4D8;** Fabiani, G., Galaris, E., Russo, L., Siettos, C., Parsimonious Physics-Informed Random Projection Neural Networks for Initial Value Problems of ODEs and index-1 DAEs. Chaos, 33, 2023. Doi: 10.1063/5.0135903 **

Last revised by G. Fabiani, March 30, 2023

We present a numerical method based on ``random projections`` with ``Gaussian kernels`` and ``physics-informed neural networks`` for the solution of initial value problems (IVP) of nonlinear STIFF ODEs and index-1 DAEs, which may also arise from the spatial discretization of PDEs.

The efficiency of the scheme is compared against two stiff ODEs/DAEs solvers, namely the ``ode15s`` and the ``ode23t`` solvers of the MATLAB ODE suite but also with DeepXDE python library (https://github.com/lululxvi/deepxde) that implements deep-learning PINNs.

Key words: ordinary differential equations (ODEs), differential algebraic equations (DAEs), stiff systems, software, Random Projection Neural Network (RPNN).

For more details and examples, see in https://doi.org/10.48550/arXiv.2203.05337

DISCLAIMER:
This software is provided "as is" without warranty of any kind.


# Matlab Examples

The main function is ada_RPNN_DAE.m (Adaptive time step solver for both ODEs and DAEs)

Here, we provide 4 examples/demos:
1) The van der Pol Equations (vdpp.m)
2) The Lotka-Volterra ODEs  (LotkaVolterra.m) [as proposed in deepxde python library]
3) The Robertson index-1 DAEs (Robertson_DAE.m)
4) A mechanics non autonomous index-1 DAEs model (ex11ode.m) [part of the benchmark problems presented in Shampine (1999), Solving Index-1 DAEs in MATLAB and Simulink]

The ODEs function handles are defined in their respective "file".m and their Jacobian in "file"_jac.m
To run the examples call "file"_main.m

# Description of the Problem
Here, we consider initial-value problems (IVPs) of ODEs and index-1 DAEs that may also arise from the spatial discretization of PDEs using for example finite differences, finite elements and spectral methods. In particular, we consider IVPs in the linear implicit form of:

$\boldsymbol{M}\dfrac{d\boldsymbol{u}(t)}{dt}  =  \boldsymbol{f}(t,\boldsymbol{u}(t)), \qquad \boldsymbol{u}(0)  =  \boldsymbol{z}. $

$\boldsymbol{u}\in \mathbb{R}^{m}$ denotes the set of the states $(u_1, u_2, \dots, u_i,\dots,  u_m)$ , $\boldsymbol{M}\in \mathbb{R}^{m\times m}$ is the so-called mass matrix with elements $M_{ij}$,
$\boldsymbol{f}: D\subseteq\mathbb{R} \times \mathbb{R}^m \to \mathbb{R}^m$ denotes a Lipschitz continuous multivariate function, with components $f_i(t, u_1,u_2,\dots,u_m)$ defined in a closed domain $D$, and $\boldsymbol{z}\in \mathbb{R}^{m}$ are the initial conditions.
When $\boldsymbol{M}=\boldsymbol{I}$, the system reduces to the canonical form.

The above formulation includes problems of DAEs when $\boldsymbol{M}$ is a singular matrix, including semi-explicit DAEs in the form:

$\dfrac{d\boldsymbol{u}(t)}{dt}  =  \boldsymbol{f}(t, \boldsymbol{u}(t),\boldsymbol{v}(t)), \quad \boldsymbol{u}(0)  =  \boldsymbol{z}, \boldsymbol{0}=\boldsymbol{g}(t,\boldsymbol{u}(t),\boldsymbol{v}(t)), $

where now $\boldsymbol{f}: \mathbb{R} \times \mathbb{R}^{m-l} \times \mathbb{R}^{l}\to \mathbb{R}^{m-l} \mbox{)}$, $\boldsymbol{g}: \mathbb{R} \times \mathbb{R}^{m-l} \times \mathbb{R}^{l}\to \mathbb{R}^{l}$ and we assume that the Jacobian $\nabla_{\boldsymbol{v}} \boldsymbol{g}$ is nonsingular.


# Documentation of the Code
  We provide an user-friendly and MATLAB-friendly software for solving ordinary differential equations (ODEs) and index-1 differential-algebraic equations (DAEs)
  using Physics-informed Random Projection Neural Network (PIRPNN). The RPNN algorithm is a fast and efficient machine learning algorithm for function approximation.
  The code supports both single and multi-dimensional systems. It allows the user to easily specify the function handle, time interval and initial condition, preserving the same structure of the MATLAB ODE suite. The user can also provide an analytical Jacobian and, in case of a DAEs, a mass matrix. 

  MATLAB CALL:
  [TOUT,YOUT] = ada_RPNN_DAE(ODEFUN,TSPAN,Y0)
 
  integrates the system of differential equations
  y' = f(t,y), from time T0 to TFINAL and y(T0)=Y0.
 
  INPUTS:
  ODEFUN: function handle
  For a row-vector of times T and a matrix Y (rows: times,colums:
  equations) ODEFUN(T,Y) must return a matrix, with columns
  corresponding to f(t,y).
  - TSPAN = [T0 TFINAL]
  To obtain solutions at specific times T0,T1,...,TFINAL (all increasing), use TSPAN =
  [T0 T1 ... TFINAL].
  - Y0: initial conditions, should include also the starting value for DAEs
 
  OUTPUTS: Each column in the solution array YOUT corresponds
  to a time returned in the row vector TOUT.
 
  [TOUT,YOUT] = ada_RPNN_DAE(ODEFUN,TSPAN,Y0,OPTIONS) solves as above with default
  integration properties replaced by values in OPTIONS (struct with fields)
  Commonly used options are relative error tolerance 'RelTol'
  (1e-3 by default) and absolute error tolerances 'AbsTol' (1e-6 by default);

  OPTIONS.Jacobian is used for giving in input the Jacobian function handle
  which for the current (and fast) implementation needs a transposed Jacobian
  e.g. J=[df_1/dx_1,df_2/dx_1; df_1/dx_2,df_2/dx_2]; that can be evaluated for
  more than 1 time istants, giving in output then a matrix, if np are the number
  of collocation points and nEq are the number of equations, of size [np x nEq*nEq].
  Otherwise the Jacobian is estimated with finite differences.

  OPTIONS.Mass is used for giving in input the Mass matrix of a DAE problem
  in the form M*dydt=f(t,y)
 
  As output one can also get a third argument for info about the computations
  [TOUT,YOUT,INFO] = ada_RPNN_DAE(ODEFUN,TSPAN,Y0,OPTIONS)
  INFO.accept: number of accepted adaptive steps
  INFO.reject: number of rejected adaptive steps
  INFO.exec_time: execution time
  INFO.num_steps: number of total step needed
  INFO.num_fun_eval: number of function evaluations
  INFO.num_Jac_solv: number of jacobian inversions(pinv/QR decompositions)

# Auxiliary functions
"Jac_Create.m" and "Jac_Create_Sparse.m" are two auxiliary subroutines, written in C (Programming Language) that can be run in MATLAB as a MEX file, for constructing the "big" Jacobian matrix (w.r.t to the weights of the RPNN) of the Newton iterations.
