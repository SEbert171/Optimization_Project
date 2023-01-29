
clear clc;

import casadi.*
syms X_sym Y_sym

X = casadi.MX(1,1);
Y = casadi.MX(1,1);

%Settings

%use full Casadi_Ipopt implementation, otherwise fall back to our
%Newton-implementation
Use_Ipopt_solver = false;

%maximum number of NLP iterations:
max_NLP_iterations = 30;

%maximum number of Newton iterations:
max_Newton_iterations = 30;

%test_function_type: ackley and rastrigin are possible
test_function_type='ackley';

%show 3dplot--> very calculation hungry but good to see whats wrong
Show3dplot=false;



%Define objective and constraint functions
% ackley function
if (strcmp(test_function_type, 'ackley'))
    f = -20*exp(-0.2*sqrt(0.5*(X^2+Y^2)))-exp(0.5*(cos(2*pi*X)+cos(2*pi*Y)))+exp(1)+20;
    f_sym = -20*exp(-0.2*sqrt(0.5*(X_sym^2+Y_sym^2)))-exp(0.5*(cos(2*pi*X_sym)+cos(2*pi*Y_sym)))+exp(1)+20;
elseif (strcmp(test_function_type, 'rastrigin'))
    f = 20+X^2-10*cos(2*pi*X)+Y^2-10*cos(2*pi*Y);
    f_sym = 20+X_sym^2-10*cos(2*pi*X_sym)+Y_sym^2-10*cos(2*pi*Y_sym);
else
    msg='Test function not recognized. Use ackley or rastrigin.';
    error(msg);
end


% constraint function
g = X^2+Y^2-25;
g_sym = X_sym^2+Y_sym^2-25;

%Casadi functions for being able to evaluate f and g
f_eval = Function('f',{[X,Y]},{f});
g_eval = Function('g',{[X,Y]},{g});


% starting point
iguess= [-2.7;-1.3]; % try -5;-5, -4;-4, -3;-3 ...
x=iguess(1);
y=iguess(2);

% starting penalty parameter
gamma=10;

% convergence criteria for newton method
e = 10^(-10);

%maximum constraint violation
max_constraint_violation = 10 ^(-4);


tic


if Use_Ipopt_solver
    [rguess, constraint_violation] = Ipopt_solve_NLPs(f,g,gamma,max_constraint_violation,max_NLP_iterations, iguess);
else
    [rguess, constraint_violation] = solve_Penalty_NLP_Newton(f,g,gamma,max_constraint_violation,max_NLP_iterations, iguess, e);
end
time_elapsed = toc;



% displaying the results
% rguess(:,n+1)=[];

% plots:
figure(1)
fcontour(f_sym, 'Fill', 'On');
hold on;
plot(rguess(1),rguess(2),'*-r');
grid on;
fimplicit(g_sym,'r');
colorbar

if(Show3dplot)
figure(2)
fmesh(f_sym)
hold on;
plot3(rguess(1,:),rguess(2,:),subs(f,[X,Y], [rguess(1),rguess(2)]),'*-r');
grid on
fimplicit3([g,0])
colorbar
end

% output:
disp(['Initial Objective Function Value: %d\n\n',full(f_eval(starting_point))]);
if (norm(Jp) < e)
 fprintf('Minimum succesfully obtained...\n\n');
end

fprintf('Number of Iterations for Convergence: %d\n\n', n);
fprintf('Point of Minima: [%f,%f]\n\n', rguess(1), rguess(2));
disp(['Objective Function Minimum Value after Optimization: %f\n\n', full(f_eval(rguess))]);

fprintf('Norm of the constraint violation: %f\n\n', norm(constraint_violation));


