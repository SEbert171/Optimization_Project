clear all
import casadi.*



%set test function
test_function_type = 'rosenbrock';
%grid distance
h = 0.1;

A(:,:,1)=zeros(10/h+1);
A(:,:,2)=zeros(10/h+1);



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




[X,Y] = meshgrid(-6:0.1:6,-6:0.1:6);
if (strcmp(test_function_type,'rosenbrock') || strcmp(test_function_type,'convex'))
    [X,Y] = meshgrid(-2:0.1:2,-2:0.1:2);
end
%Definition of the different functions as matlab functions for plotting
%(equivalent with the definition in main
if (strcmp(test_function_type, 'ackley'))
    f1 = -20*exp(-0.2*sqrt(0.5*(X.^2+Y.^2 +10^(-3))))-exp(0.5*(cos(2*pi*X)+cos(2*pi*Y)))+exp(1)+20;
    g1 = X.^2+Y.^2-25;% constraint function for ackley
elseif (strcmp(test_function_type, 'rastrigin'))
    f1 = 20+X.^2-10*cos(2*pi*X)+Y.^2-10*cos(2*pi*Y); % rastrigin function
    g1 = X.^2+Y.^2-26.2144; % constraint function for rastrigin
elseif (strcmp(test_function_type, 'rosenbrock'))
    f1 = (1-X).^2+100*(Y-X.^2).^2;
    g1 = X.^2+Y.^2-1.5; % constraint function for Rosenbrock
elseif (strcmp(test_function_type, 'convex'))
    f1 = X.^2+X.*Y+Y.^2+exp(X);
    g1 = X.^2+Y.^2-1;
else
    msg='Test function not recognized. Use ackley, rastrigin, rosenbrock or convex.';
    error(msg);
end


% plots:
figure(8)
contour(X,Y,f1, 'Fill', 'On');
hold on

for i = 1:4/h+1
    for j = 1:4/h+1
        x = h*i-2-h;
        y = h*j-2-h;
        [solution,~, gamma] = main_test([x;y], test_function_type, 10^(-4), (10^-6));
        A(i,j,:) = solution;
        if (norm(full(g(solution(1),solution(2)))) <= 10^(-4) && check_SOSC(solution,test_function_type,gamma,10^(-4)))
            plot(x,y,'*g')%plot green if it did converge to the constraint and the SOSC is satisfied
            hold on
            converged = true;
        elseif (norm(full(g(solution(1),solution(2)))) <= 10^(-4))
            plot(x,y,'*y')%plot yellow if it did converge to the constraint
            hold on
            converged = false;
        else
            plot(x,y,'*r')%plot red if it did not converge to the constraint
            hold on
            converged = false;
        end
        disp(['SOLVED PROBLEM FOR [X,Y]=[',num2str(x),',',num2str(y),']'])
    end
end

time_Newton = ones(1,7);
time_Ipopt = ones(1,7);
solution_Newton = ones(2,7);
solution_Ipopt = ones(2,7);
difference = ones(1,7);
for i = 2:8
    [solution_Newton(:,i-1),time_Newton(i-1)] = main_test([2;3],'rastrigin',10^(-i), 10^(-6));
    [solution_Ipopt(:,i-1),time_Ipopt(i-1)] = maini([2;3],'rastrigin',10^(-i));
end

figure(2)
semilogy(2:8,time_Newton,'-g')
hold on
semilogy(2:8,time_Ipopt,'-r')
legend('Quadratic Penalty method','IpOpt')
title('Computation time comparison')
xlabel('exponent of max constraint violation')
ylabel('Required Computation time in seconds')




