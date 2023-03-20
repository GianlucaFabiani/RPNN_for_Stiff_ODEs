# RPNN_for_Stiff_ODEs
Physics Informed Random Projection Neural Network for solving Stiff ODEs and DAEs

We present a numerical method based on ``random projections`` with ``Gaussian kernels`` and ``physics-informed neural networks`` for the solution of initial value problems (IVP) of nonlinear ODEs and index-1 DAEs, which may also arise from the spatial discretization of PDEs. The internal weights are fixed to ones while the unknown weights between the hidden and output layer are computed with Newton's iterations, using the Moore-Penrose pseudo-inverse for low to medium scale, and sparse QR decomposition with $L^2$ regularization for medium to large scale systems. Building on previous works on random projections, we also prove its approximation accuracy. To deal with stiffness and sharp gradients, we propose an adaptive step-size scheme, and address a continuation method for providing good initial guesses for the Newton iterations. The "optimal" bounds of the uniform distribution from which the values of the shape parameters of the Gaussian kernels are sampled and the number of basis functions are "parsimoniously" chosen based on a bias-variance trade-off decomposition. To assess the performance of the scheme in terms of both numerical approximation accuracy and computational cost, we used eight benchmark problems (three index-1 DAEs problems, and five stiff ODEs problems including the Hindmarsh-Rose neuronal model of chaotic dynamics and the Allen-Cahn phase-field PDE). The efficiency of the scheme was compared against two stiff ODEs/DAEs solvers, namely the ``ode15s`` and the ``ode23t`` solvers of the MATLAB ODE suite as well as against deep learning as implemented in the ``DeepXDE`` library for scientific machine learning and physics-informed learning, for the solution of the Lotka-Volterra ODEs included in the demos of the library.

https://doi.org/10.48550/arXiv.2203.05337
