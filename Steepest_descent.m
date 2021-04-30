function [xstar, fval,z,a,f_cont,i] = Steepest_descent(f,x,tol,iter,method)
% Steepest descent algorithm
fprintf('Minimization using Steepest Descent algorithm: \n');
x0 = x;
z = x';
alpha0 = 0.002;          % Steepest decsent     
c1 = 0.0001; 
c2 = 0.1;                % Steepest decsent  
a = alpha0;
f_cont = 0;
fval = [];
fprintf('%12s %12s %8s %16s\n','Iteration','Func-cont','f(x)','Step-size');
for i = 1:iter
    g = grad(x);
    if method == 1
        [alpha, f_iter] = Backtracking_Armijo(f,g,x,-g,alpha0);
    elseif method == 2
        [alpha, f_iter] = Wolfe_condition(f,x,-g,alpha0,c1,c2);
    else
        [alpha, f_iter] = Strong_Wolfe(f,x,-g,alpha0,c1,c2);
    end
    f_cont = f_cont+f_iter;
    a = [a,alpha];
    xstar = x - alpha*g;
    z = [z;xstar'];
    fval = [fval,f(xstar)];
    gnorm = norm(grad(xstar));
    if ~isfinite(xstar)
        fprintf(1,'Number of iterations: %d\n', i);
        error('x is inf or NaN')
    end
    if gnorm < tol*max(1,norm(grad(x0)))
        fprintf(1,'Convergance criterion: %d\n', tol);
        fprintf(1,'Maximum gradient component: %d\n', gnorm);
        fprintf(1,'Number of iterations: %d\n', i-1);
        break;
    end
    x = xstar;
    fprintf('%10d %10d %14.2e %14.2e\n',i,f_iter,fval(i),alpha);
end
if i == iter
    fprintf(1,'Number of iterations is: %d\n', iter);
    fprintf(1,'Convergance criterion: %d\n', tol);
    fprintf(1,'Maximum gradient component: %d\n', gnorm);
end