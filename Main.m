clc; clear;
%%  Initial variables
x = [-1.2; 1];

% select a linesearch methods:
method = 3; %      1: Backtracking_Armijo      
            %      2: Wolfe condition
            %      3: Strong Wolfe conditions

%% optimization
tol = 10^(-8);
iter = 20000;

% Active one of the following linesearch direction algorithms:
%%%%%% Steepest descent
% [Epar,fval,z,alpha,f_cont,iter] = Steepest_descent(@Rosenbrock,x,tol,iter,method);
%%%%%% Modified Newton
[Epar,fval,z,alpha,f_cont,iter] = Modified_Newton(@Rosenbrock,x,tol,iter,method);
% %%%%% BFGS
%[Epar,fval,z,alpha,f_cont,iter] = BFGS(@Rosenbrock,x,tol,iter,method);

%% result
fprintf(1,'Function evaluation: %d\n', f_cont);
fprintf(1,'fval: %e\n', fval(end));
fprintf(1,'Estimated x1 is: %d\n', Epar(1));
fprintf(1,'Estimated x2 is: %d\n', Epar(2));

%% Illustration
syms x1 x2
f = 100*(x2 - x1^2 )^2 + (1 - x1)^2;
fcontour(f,[-1.2 2.2 -1.5 1.5])
hold on
figure(1), plot(z(:,1),z(:,2),'r*-') 
xlabel('x1'); ylabel('x2'); 
figure(2), plot(1:length(alpha),alpha,'LineWidth',1.5)
xlabel('Iteration');  ylabel('alpha'); ylim([0,1.2]); box off
figure(3), plot(1:length(fval),log(abs(fval)),'LineWidth',1.5)
xlabel('Iteration');  ylabel('log(f(x))'); box off

%% fminunc
fprintf('Minimization MATLAB fminunc with Quasi-Newton algorithm: \n');
history= histclass;
outf= @(x,optimValues,state)outfun(x,optimValues,state,history);
options = optimoptions(@fminunc,'Display','iter','Algorithm',...
         'quasi-newton','HessUpdate','bfgs','OutputFcn',outf,...
         'TolFun',1e-8,'MaxFunEvals',2e4,'MaxIter',2e4);
xfmin = fminunc(@Rosenbrock,x,options);
fprintf('Number of function evaluations: %d\n', history.count);
hold on
plot(1:length(history.fval),log(history.fval),'LineWidth',1.5)

%% ADMB
dat = 'ADMB.txt';  % fun values
z = load(dat);
plot(1:length(z),log(z),'LineWidth',1.5)
legend('My result','fminunc','ADMB')

%%
%close all;