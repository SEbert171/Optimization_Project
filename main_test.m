function [solution, time_elapsed, gamma] = main_test(starting_point, test_function_type, max_constraint_violation, Newton_terminal_condition)
%Copy of the test function to implement multiple starting point analysis

import casadi.*
import casadi.*

%Settings
%maximum number of NLP iterations:
max_NLP_iterations =50;
%maximum number of Newton iterations:
max_Newton_iterations = 100;
max_line_search_iterations = 10;
%search method type: exact_Newton or Ipopt are possible
solve_method = 'exact_Newton';


%Casadi initialization
X = MX.sym('X',1);
Y = MX.sym('Y',1);

if (strcmp(test_function_type, 'ackley'))
    f = Function('ackley',{X,Y},{-20*exp(-0.2*sqrt(0.5*(X^2+Y^2 +10^(-3))))-exp(0.5*(cos(2*pi*X)+cos(2*pi*Y)))+exp(1)+20});
    g = Function('constraint',{X,Y},{X^2+Y^2-25}); % constraint function for ackley
elseif (strcmp(test_function_type, 'rastrigin'))
    f = Function('rastrigin',{X,Y},{20+X^2-10*cos(2*pi*X)+Y^2-10*cos(2*pi*Y)}); % rastrigin function
    g = Function('constraint',{X,Y},{X^2+Y^2-26.2144}); % constraint function for rastrigin
elseif (strcmp(test_function_type, 'rosenbrock'))
    f = Function('rosenbrock',{X,Y},{(1-X)^2+100*(Y-X^2)^2});
    g = Function('constraint',{X,Y},{X^2+Y^2-1.5}); % constraint function for Rosenbrock
elseif strcmp(test_function_type, 'convex')
    f = Function('convex',{X,Y},{X^2+Y^2+X*Y+exp(X)});
    g = Function('constraint',{X,Y},{X^2+Y^2-1});
else
    msg='Test function not recognized. Use ackley, rastrigin, rosenbrock or convex.';
    error(msg);
end



iguess = starting_point;
solution_points = iguess;
constraint_violation = norm(full(evalf(g(iguess(1),iguess(2)))));
% penalty parameter
gamma=10;

rguess=iguess; %rguess -> record guesses; iguess -> initial guess

n=1;

tic
F = Function('Penalty',{X,Y},{f(X,Y)+0.5*gamma*g(X,Y).^2});
Q = Function('Penalty_part',{X,Y},{0.5*gamma*g(X,Y).^2});
[Jp, Hp] = calculate_derivatives(F,rguess);
while (norm(constraint_violation(length(constraint_violation))) > max_constraint_violation && n<=max_NLP_iterations && norm(Jp) > Newton_terminal_condition*10^-n)
    disp(['NLP-iteration: ',num2str(n)]);
        if strcmp(solve_method,'exact_Newton')
            rguess = solve_Penalty_NLP_Newton(F,Q,iguess,Newton_terminal_condition,max_Newton_iterations, max_line_search_iterations);
        elseif strcmp(solve_method, 'Ipopt')%Solves the NLP directly with Ipopt
            rguess = solve_Penalty_NLP_IpOpt(test_function_type, gamma, iguess);
        else
            error('Search direction method not recognized. Use exact_Newton, constraint_Newton or Ipopt.')
        end
        
    solution_points = [solution_points rguess];
    constraint_violation = [constraint_violation, norm(full(g(rguess(1),rguess(2))))];
        if (norm(constraint_violation(length(constraint_violation))) <= 10^(-20-n+1))
            break
        end

    n=n+1;
    gamma=gamma*1.5;
    F = Function('Penalty',{X,Y},{f(X,Y)+0.5*gamma*g(X,Y)^2});
    iguess = rguess;
    [Jp, Hp] = calculate_derivatives(F,iguess);
    
end
time_elapsed = toc;


solution = rguess;
end