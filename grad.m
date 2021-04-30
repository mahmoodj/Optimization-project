function g = grad(x)
%  Rosenbrock gradient

g = [2*x(1) - 400*x(1)*(- x(1)^2 + x(2)) - 2; - 200*x(1)^2 + 200*x(2)];