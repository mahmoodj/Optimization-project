function [alpha, l] = Wolfe_condition(f,x,p,alpha,c1,c2)
maxiter = 1e2;
tau = 0.5;
l = 1;
while 1
    if f(x+alpha*p) <= f(x)+c1*alpha*grad(x)'*p && grad(x+alpha*p)'*p >= c2*grad(x)'*p
        break;
    end
    alpha = tau*alpha;
    l = l+1;
    if l>maxiter
        alpha = 0.0005;
        break;
    end
end
if l~=1
    l = l-1;
   alpha = alpha/tau;
end