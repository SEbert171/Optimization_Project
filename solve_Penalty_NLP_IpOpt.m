function [rguess] = solve_Penalty_NLP_IpOpt(test_function_type, gamma, iguess)
%Solves the Penalty NLP using the Ipopt solver. F is the penalty function (f+0.5*gamma*g^2) and
%iguess is the initial guess for the Ipopt solver. Returns the new
%approximation for a minimum
import casadi.*
opti = casadi.Opti();

X_ipopt = opti.variable();
Y_ipopt = opti.variable();

if (strcmp(test_function_type, 'ackley'))
    f_ipopt = -20*exp(-0.2*sqrt(0.5*(X_ipopt^2+Y_ipopt^2+1e-6)))-exp(0.5*(cos(2*pi*X_ipopt)+cos(2*pi*Y_ipopt)))+exp(1)+20; % ackley function
    g_ipopt = X_ipopt^2+Y_ipopt^2-25; % constraint function for ackley
elseif (strcmp(test_function_type, 'rastrigin'))
    f_ipopt = 20+X_ipopt^2-10*cos(2*pi*X_ipopt)+Y_ipopt^2-10*cos(2*pi*Y_ipopt); % rastrigin function
    g_ipopt = X_ipopt^2+Y_ipopt^2-26.2144; % constraint function for rastrigin
elseif (strcmp(test_function_type, 'rosenbrock'))
    f_ipopt = (1-X_ipopt)^2+100*(Y_ipopt-X_ipopt^2)^2;
    g_ipopt = X_ipopt^2+Y_ipopt^2-1.5; % constraint function for Rosenbrock 
end


opti.minimize(f_ipopt+0.5*gamma*g_ipopt^2)
opti.solver('ipopt');
opti.set_initial(X_ipopt, iguess(1));
opti.set_initial(Y_ipopt, iguess(2));
sol = opti.solve();
 
 x=sol.value(X_ipopt);
 y=sol.value(Y_ipopt);
 rguess=[x;y];

end