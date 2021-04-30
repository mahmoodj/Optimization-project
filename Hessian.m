function h = Hessian(x)
%  Rosenbrock hessian

h = 200*[6*x(1)^2 - 2*x(2) + 0.01, - 2*x(1) ; - 2*x(1), 1];