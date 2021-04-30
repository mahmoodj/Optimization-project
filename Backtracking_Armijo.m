function [alpha, l] = Backtracking_Armijo(f,g,x,p,alpha)
maxiter = 1e2;
tau = 0.5; 
eta = 0.1;
l = 1;
while f(x+alpha*p) > f(x)+eta*alpha*g'*p
    alpha = tau*alpha;
    l = l+1;
    if l>maxiter
        alpha = 0.002;
        break;
    end
end
if l~=1
    l = l-1;
   alpha = alpha/tau;
end