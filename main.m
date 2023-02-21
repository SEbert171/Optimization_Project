clear all;
import casadi.*


%Settings
%maximum number of NLP iterations:
max_NLP_iterations =50;
%maximum number of Newton iterations:
max_Newton_iterations = 100;
max_line_search_iterations = 10;
%test_function_type: ackley and rastrigin are possible
test_function_type='ackley';
%search method type: exact_Newton or constraint_Newton are possible
solve_method = 'exact_Newton';
%show 3dplot--> very calculation hungry but good to see whats wrong
Show3dplot=false;


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
else
    msg='Test function not recognized. Use ackley, rastrigin or rosenbrock.';
    error(msg);
end

% starting point
starting_point= [-1.8;-2]; % try -5;-5, -4;-4, -3;-3 ...
iguess = starting_point;
x=iguess(1);
y=iguess(2);
solution_points = iguess;
constraint_violation = norm(full(evalf(g(iguess(1),iguess(2)))));
% penalty parameter
gamma=10;
% convergence criteria for newton method
Newton_terminal_condition = 10^(-6);
% %maximum constraint violation
max_constraint_violation = 10 ^(-4);
% tic

rguess=iguess; %rguess -> record guesses; iguess -> initial guess

n=1;

tic
F = Function('Penalty',{X,Y},{f(X,Y)+0.5*gamma*g(X,Y)^2});
Q = Function('Penalty_part',{X,Y},{0.5*gamma*g(X,Y)^2});
[Jp, Hp] = calculate_derivatives(F,rguess);
while (norm(constraint_violation(length(constraint_violation))) > max_constraint_violation && n<=max_NLP_iterations && norm(Jp) > Newton_terminal_condition*10^-n)
    disp(['NLP-iteration: ',num2str(n)]);
        if strcmp(solve_method,'exact_Newton')
            rguess = solve_Penalty_NLP_Newton(F,Q,iguess,Newton_terminal_condition,max_Newton_iterations, max_line_search_iterations);
        elseif strcmp(solve_method,'constraint_Newton')% really include?
            error('Not yet implemented: Change calculation in the function "solve_Penalty_NLP_augmented_Newton"')
            rguess = solve_Penalty_NLP_augmented_Newton(F,g,Q,gamma,iguess,Newton_terminal_condition,max_Newton_iterations,max_line_search_iterations);
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


    %Try of including SOSC
    %[Jp,Hp] = penalty_derivatives(F,iguess);
    %for i = 1:length(Hp(1,:))
         %if (norm(Jp)<10^-4 && eigenvalues_Hp(i)<0) %Backstop if SOSC is not fulfilled, at the moment via change of gamma, but better find new iguess
              %iguess = iguess+[0.2;0.2];
              %gamma = gamma/1000;
              %disp(['SOSC not fulfilled. Set gamma to ',num2str(gamma)])
              %break
         %end
    %end
    n=n+1;
    gamma=gamma*10;
    F = Function('Penalty',{X,Y},{f(X,Y)+0.5*gamma*g(X,Y)^2});
    iguess = rguess;
    [Jp, Hp] = calculate_derivatives(F,iguess);
end
time_elapsed = toc;


%plotting the results
plots(solution_points, constraint_violation,test_function_type)


% output of the results:
fprintf('Initial Objective Function Value: %d\n\n',full(f(starting_point(1),starting_point(2))));
if (norm(Jp) < Newton_terminal_condition)
 fprintf('Minimum succesfully obtained...\n\n');
end

fprintf('Number of Iterations for Convergence: %d\n\n', n);
fprintf('Point of Minima: [%f,%f]\n\n', rguess(1), rguess(2));
fprintf('Objective Function Minimum Value after Optimization: %f\n\n', full(f(rguess(1),rguess(2))));

fprintf('Norm of the constraint violation: %f\n\n', norm(constraint_violation(length(constraint_violation))));
disp('Solution points were:')
disp(solution_points)
disp(['Solution time was ',num2str(time_elapsed),' seconds'])