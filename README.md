# RPNN_for_Stiff_ODEs
Physics Informed Random Projection Neural Network for solving Stiff ODEs and DAEs

We present a numerical method based on ``random projections`` with ``Gaussian kernels`` and ``physics-informed neural networks`` for the solution of initial value problems (IVP) of nonlinear ODEs and index-1 DAEs, which may also arise from the spatial discretization of PDEs.

The efficiency of the scheme is compared against two stiff ODEs/DAEs solvers, namely the ``ode15s`` and the ``ode23t`` solvers of the MATLAB ODE suite.

https://doi.org/10.48550/arXiv.2203.05337


# Maltlab Code

The main function is ada_RPNN_DAE (Adaptive time step solver for both ODEs and DAEs)

Here we have 4 example:
1) van der Pol Equations
2) Lotka-Volterra ODEs
3) Robertson DAEs
4) 

The ODEs function handles are defined in their respective <file>.m and their Jacobian in <file>_jac.m
To run the examples call <file>_main.m
