function [xstar, fval,z,a,f_cont,i] = BFGS(f,x,tol,iter,method)
% BFGS Quasi-Newton algorithm
fprintf('Minimization using BFGS Quasi-Newton algorithm: \n');
m = 1;
x0 = x;
z = x';
alpha = 0.007;
a = alpha;
f_cont = 0;
fval = f(x);
for i = 1:m
    g = grad(x);
    gnorm = norm(g);
    xstar = x - alpha*g;
    z = [z;xstar'];
    %fval = [fval, f(xstar)];
    if ~isfinite(xstar)
        fprintf(1,'Number of iterations: %d\n', i);
        error('x is inf or NaN')
    end
    if gnorm < tol*max(1,norm(grad(x0)))
        fprintf(1,'Gradient is less than: %d\n', tol);
        fprintf(1,'Number of iterations: %d\n', i);
        break;
    end
    xold = x;
    x = xstar;
end
fprintf('%12s %12s %8s %16s\n','Iteration','Func-cont','f(x)','Step-size');
fprintf('%10d %10d %14.2e %14.2e\n',1,1,fval(i),alpha);
alpha0 = 0.9;           % BFGS
c1 = 0.0001; 
c2 = 0.9;               % Newton, BFGS
xnew = xstar;
nvar = length(x);
B = eye(nvar);
g = grad(xnew);
for i = 1:iter
    y = grad(xnew) - grad(xold);
    s = xnew - xold;
    B = B - (1/(s'*B*s))*(B*s)*(B*s)' + (1/(y'*s))*(y*y');
    p = -inv(B)*g;
    if method == 1
        [alpha, f_iter] = Backtracking_Armijo(f,g,xnew,p,alpha0);
    elseif method == 2
        [alpha, f_iter] = Wolfe_condition(f,xnew,p,alpha0,c1,c2);
    else
        [alpha, f_iter] = Strong_Wolfe(f,xnew,p,alpha0,c1,c2);
    end
    f_cont = f_cont+f_iter;
    a = [a,alpha];
    xstar = xnew + alpha*p;
    z = [z;xstar'];
    g = grad(xstar);
    gnorm = norm(g);
    fval = [fval,f(xstar)];
    fprintf('%10d %10d %14.2e %14.2e\n',i+m,f_iter,f(xstar),alpha);
    if ~isfinite(xstar)
        fprintf(1,'Number of iterations: %d\n', i);
        error('x is inf or NaN')
    end
    if gnorm < tol*max(1,norm(grad(x0)))
        fprintf(1,'Convergance criterion: %d\n', tol);
        fprintf(1,'Maximum gradient component: %d\n', gnorm);
        fprintf(1,'Number of iterations: %d\n', i);
        break;
    end
    xold = xnew;
    xnew = xstar;
end
if i == iter
    fprintf(1,'Number of iterations is: %d\n', iter);
    fprintf(1,'Convergance criterion: %d\n', tol);
    fprintf(1,'Maximum gradient component: %d\n', gnorm);
end