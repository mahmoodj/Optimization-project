function [xstar, fval,z,a,f_cont,i] = Modified_Newton(f,x,tol,iter,method)
% Modified Newton algorithm
fprintf('Minimization using Modified Newton algorithm: \n');
x0 = x;
nvar = length(x);
beta = 1e4;
z = x';
alpha0 = 1.0;           % Newton 
c1 = 0.0001; 
c2 = 0.9;               % Newton, BFGS
a = alpha0;
f_cont = 0;
fval = [];
fprintf('%12s %12s %8s %16s\n','Iteration','Func-cont','f(x)','Step-size');
for i = 1:iter
    B = Hessian(x);
    [V,D] = eig(B);
    if nnz(B)== 0
        epsilon = 1;
    else
        epsilon = norm(B)/beta;
    end
    for j=1:nvar
        if D(j,j) >= epsilon
            D(j,j) = D(j,j);
        elseif D(j,j) <= -epsilon
            D(j,j) = -D(j,j);
        else
            D(j,j) = epsilon;
        end
    end
    B = V*D*V';
    g = grad(x);
    p = -inv(B)*g;
    if method == 1
        [alpha, f_iter] = Backtracking_Armijo(f,g,x,p,alpha0);
    elseif method == 2
        [alpha, f_iter] = Wolfe_condition(f,x,p,alpha0,c1,c2);
    else
        [alpha, f_iter] = Strong_Wolfe(f,x,p,alpha0,c1,c2);
    end
    f_cont = f_cont+f_iter;
    a = [a,alpha];
    xstar = x + alpha*p;
    z = [z;xstar'];
    g = grad(xstar);
    gnorm = norm(g);
    fval = [fval,f(xstar)];
    fprintf('%10d %10d %14.2e %14.2e\n',i,f_iter,fval(i),alpha);
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
    x = xstar;
end
if i == iter
    fprintf(1,'Number of iterations is: %d\n', iter);
    fprintf(1,'Convergance criterion: %d\n', tol);
    fprintf(1,'Maximum gradient component: %d\n', gnorm);
end