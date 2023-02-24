function[solution, time_elapsed] = maini(starting_point, test_function_type, max_constraint_violation)

import casadi.*



opti = casadi.Opti();

X = opti.variable();
Y = opti.variable();

%Definition of the different functions as matlab functions for plotting
%(equivalent with the definition in main
if (strcmp(test_function_type, 'ackley'))
    f = -20*exp(-0.2*sqrt(0.5*(X.^2+Y.^2 +10^(-3))))-exp(0.5*(cos(2*pi*X)+cos(2*pi*Y)))+exp(1)+20;
    g = X.^2+Y.^2-25;% constraint function for ackley
elseif (strcmp(test_function_type, 'rastrigin'))
    f = 20+X.^2-10*cos(2*pi*X)+Y.^2-10*cos(2*pi*Y); % rastrigin function
    g = X.^2+Y.^2-26.2144; % constraint function for rastrigin
elseif (strcmp(test_function_type, 'rosenbrock'))
    f = (1-X).^2+100*(Y-X.^2).^2;
    g = X.^2+Y.^2-1.5; % constraint function for Rosenbrock
elseif (strcmp(test_function_type, 'convex'))
    f = X.^2+X.*Y+Y.^2+exp(X);
    g = X.^2+Y.^2-1;
else
    msg='Test function not recognized. Use ackley, rastrigin, rosenbrock or convex.';
    error(msg);
end

tic
opti.minimize(f)
 
opti.set_initial(X, starting_point(1));
opti.set_initial(Y, starting_point(2));
%opti.constr_viol_tol = max_constraint_violation;

opti.subject_to(g==0)

p_opts = struct();
		s_opts = struct("constr_viol_tol", max_constraint_violation);

opti.solver('ipopt', p_opts, s_opts);
sol = opti.solve();
 
x=sol.value(X);
y=sol.value(Y);
time_elapsed = toc;

disp('Solution point is:')
disp([x;y]);

disp(['Solution time was ',num2str(time_elapsed),' seconds'])
solution = [x;y];

