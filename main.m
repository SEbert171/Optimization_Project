clear all;
syms X Y
import casadi.*


%Settings
%maximum number of NLP iterations:
max_NLP_iterations =7;
%maximum number of Newton iterations:
max_Newton_iterations = 10;
max_line_search_iterations = 5;
%test_function_type: ackley and rastrigin are possible
test_function_type='ackley';
%search method type: exact_Newton or constraint_Newton are possible
solve_method = 'exact_Newton';
%show 3dplot--> very calculation hungry but good to see whats wrong
Show3dplot=false;

if (strcmp(test_function_type, 'ackley'))
    f = -20*exp(-0.2*sqrt(0.5*(X^2+Y^2+1e-6)))-exp(0.5*(cos(2*pi*X)+cos(2*pi*Y)))+exp(1)+20; % ackley function
    g = X^2+Y^2-25; % constraint function for ackley
elseif (strcmp(test_function_type, 'rastrigin'))
    f = 20+X^2-10*cos(2*pi*X)+Y^2-10*cos(2*pi*Y); % rastrigin function
    g = X^2+Y^2-26.2144; % constraint function for rastrigin
elseif (strcmp(test_function_type, 'rosenbrock'))
    f = (1-X)^2+100*(Y-X^2)^2;
    g = X^2+Y^2-1.5; % constraint function for Rosenbrock
else
    msg='Test function not recognized. Use ackley or rastrigin.';
    error(msg);
end

% starting point
starting_point= [1;1]; % try -5;-5, -4;-4, -3;-3 ...
iguess = starting_point;
x=iguess(1);
y=iguess(2);
solution_points = iguess;
constraint_violation = norm(subs(g,[X,Y],[iguess(1),iguess(2)]));
% penalty parameter
gamma=10;



% convergence criteria for newton method
e = 10^(-6);

% %maximum constraint violation
% max_constraint_violation = 10 ^(-4);
% tic
F=f+0.5*gamma*(g^2);
% jacobian and hessian derivatives of penalty function F=f+0.5*gamma*(g^2)
[Jp,Hp] = penalty_derivatives(F,iguess);

rguess=iguess; %rguess -> record guesses; iguess -> initial guess

n=0;
k=1;

tic

while (norm(constraint_violation(length(constraint_violation)))>10^(-6) && n<=max_NLP_iterations)
    l=0;
    while (double(norm(penalty_derivatives(F,iguess))) > e && l<= max_Newton_iterations)

        % backtracking with armijo condition
        t=1;
        kk = 0;
        
        if strcmp(solve_method,'exact_Newton')
            disp(['Newton iteration: ',num2str(l),' for NLP-iteration ',num2str(n)]);
            [Jp,Hp] = penalty_derivatives(F,iguess);
            search_direction = (-Hp\Jp);
                while (subs(F,[X;Y],(iguess+t*search_direction))>=double(subs(F,[X;Y],(iguess)+0.1*(Jp.')*search_direction))  && kk<=max_line_search_iterations )
                t=0.8*t;
                kk = kk+1;
                end
                disp(['Step length was: ',num2str(t)])
            iguess = iguess + t*search_direction;
            %eigenvalues_Hp = eig(Hp); %Tried to implement SOSC did not work yet
            %for i = 1:length(Hp(1,:))
                %if (norm(Jp)<10^-4 && eigenvalues_Hp(i)<0) %Backstop if SOSC is not fulfilled, at the moment via change of gamma, but better find new iguess
                    %iguess = starting_point+[0.5;-0.5];
                    %gamma = gamma/10000;
                    %disp(['SOSC not fulfilled. Set gamma to ',num2str(gamma)])
                    %break
               %else
                    rguess=iguess;
                %end
            %end
            l=l+1;
        elseif strcmp(solve_method,'constraint_Newton')% really include?
            disp(['Newton iteration: ',num2str(l),' for NLP-iteration ',num2str(n)]);
            [Jp,Hp] = penalty_derivatives(F,iguess);
            [Jg,Hg] = penalty_derivatives(g,iguess);
            [JQ, ~] = penalty_derivatives(0.5*gamma*g^2,iguess);
            B = [Hp+1/gamma*subs(g,[X,Y],[iguess(1),iguess(2)])*Hg Jg;Jg.' -gamma];
            c = [-JQ;0];
            solve_linear_system = B\c;
            search_direction = solve_linear_system(1:2);
            while (subs(F,[X;Y],(iguess+t*search_direction))>=double(subs(F,[X;Y],(iguess)+0.1*(Jp.')*search_direction))  && kk<=max_line_search_iterations )
                t=0.8*t;
                kk = kk+1;
            end
            disp(['Step length was: ',num2str(t)])
            iguess = iguess + t*search_direction;
            rguess=iguess;
            l=l+1;
        elseif strcmp(solve_method, 'Ipopt')
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
            F_ipopt = f_ipopt + 0.5*gamma*g_ipopt^2;
            opti.minimize(F_ipopt)
 
            opti.set_initial(X_ipopt, iguess(1));
            opti.set_initial(Y_ipopt, iguess(2));

            opti.solver('ipopt');
            sol = opti.solve();
 
            x=sol.value(X_ipopt);
            y=sol.value(Y_ipopt);

            disp([x;y]);
            rguess = [x;y];
            break
        else
            error('Search direction method not recognized. Use exact_Newton, constraint_Newton or Ipopt.')
        end
        
    solution_points = [solution_points rguess];
    constraint_violation = [constraint_violation, norm(subs(g,[X,Y],[rguess(1),rguess(2)]))];
        if (norm(constraint_violation) <= 10^(-20-k+1))
            break
        end

    end
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
    k=k+1;
    gamma=gamma*10;
    F=f+0.5*gamma*(g^2);
    [Jp,Hp] = penalty_derivatives(F,rguess);
end
time_elapsed = toc;
% displaying the results
% rguess(:,n+1)=[];

% plots:
figure(1)
fcontour(f, 'Fill', 'On');
hold on;
plot(solution_points(1,:),solution_points(2,:),'*-y');
grid on;
fimplicit(g,'r');
colorbar
title('first convergence steps of the algorithm')




 %plot with meshgrid
 if Show3dplot
    [X1,Y1] = meshgrid(-5:0.1:5,-5:0.1:5);
    Z = -20*exp(-0.2*sqrt(0.5*(X1^2+Y1^2+1e-6)))-exp(0.5*(cos(2*pi*X1)+cos(2*pi*Y1)))+exp(1)+20; % ackley function
    g1 = (X1.^2+Y1.^2-25).^2; % constraint function for ackley
    %figure(2)
    %surf(X1,Y1,Z)
    %hold on
    %fimplicit3([g,0])
    %hold off
    figure(3)
    surf(X1,Y1,g1)%plot of g
    title('Plot of the penalty function 0.5*mu*|g|^2')
    figure(4)
    subplot(2,2,1)
    gamma = 0.01;
    Z1 = Z+gamma*g1;
    surf(X1,Y1,Z+Z1)%plot of F
    title('mu=0.01')
    subplot(2,2,2)
    gamma = 0.1;
    Z1 = Z+gamma*g1;
    surf(X1,Y1,Z+Z1)%plot of F
    title('mu=0.1')
    subplot(2,2,3)
    gamma = 1;
    Z1 = Z+gamma*g1;
    surf(X1,Y1,Z+Z1)%plot of F
    title('mu=1')
    subplot(2,2,4)
    gamma = 10;
    Z1 = Z+gamma*g1;
    surf(X1,Y1,Z+Z1)%plot of F
    title('mu=10')
    figure(6)
    surf(X1,Y1,Z)
    title('Plot of the objective function f')
 end





 figure(5)
 semilogy(0:length(constraint_violation)-1,constraint_violation,"-o")
 title('Convergence of the constraint violation')
 xlabel('NLP-steps')
 ylabel('constraint violation')

 %Only for solution point [0,5]!!
 %error = ones(1,length(solution_points));
 %difference = ones(2,length(solution_points));
%for i = 1:length(solution_points)
%difference(:,i) = solution_points(:,i)-[0;5];
%error(i)= norm([difference(1,i) difference(2,i)]);
%end
 %figure(7)
 %semilogy(0:length(solution_points)-1,error,"-o")
 %title('Convergence of the error from the solution')
 %xlabel('NLP-steps')
 %ylabel('error')

% output:
fprintf('Initial Objective Function Value: %d\n\n',subs(f,[X,Y], [x,y]));
if (norm(Jp) < e)
 fprintf('Minimum succesfully obtained...\n\n');
end

fprintf('Number of Iterations for Convergence: %d\n\n', n);
fprintf('Point of Minima: [%f,%f]\n\n', iguess(1), iguess(2));
fprintf('Objective Function Minimum Value after Optimization: %f\n\n', subs(f,[X,Y], [iguess(1),iguess(2)]));

fprintf('Norm of the constraint violation: %f\n\n', norm(constraint_violation(length(constraint_violation))));
disp('Solution points were:')
disp(solution_points)
disp(['Solution time was ',num2str(time_elapsed),' seconds'])