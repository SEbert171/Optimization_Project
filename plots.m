function [] = plots(solution_points,constraint_violation,test_function_type, gamma_real)
%Plots the results when the algortihm has found a solution

%Settings
Show3dplot = false;


[X,Y] = meshgrid(-6:0.1:6,-6:0.1:6);
if strcmp(test_function_type,'rosenbrock')
    [X,Y] = meshgrid(-2:0.1:2,-2:0.1:2);
end
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
else
    msg='Test function not recognized. Use ackley, rastrigin or rosenbrock.';
    error(msg);
end



% plots:
figure(1)
contour(X,Y,f, 'Fill', 'On');
hold on;
plot(solution_points(1,:),solution_points(2,:),'*-y');
grid on;
theta = linspace(0,2*pi);
x = 5*cos(theta);
y = 5*sin(theta);
if strcmp(test_function_type,'rosenbrock')
    x = sqrt(1.5)*cos(theta);
    y = sqrt(1.5)*sin(theta);
end
plot(x,y,'-r')
colorbar
title('first convergence steps of the algorithm')




 %plot with meshgrid
 if Show3dplot
    figure(3)
    surf(X,Y,g)%plot of g
    title('Plot of g')
    figure(4)
    subplot(2,2,1)
    gamma = 0.01;
    F = f+0.5*gamma*g.^2;%.^2;
    surf(X,Y,F)%plot of F
    title('Plot of F=f+0.5*mu*g^2 for mu=0.01')
    subplot(2,2,2)
    gamma = 0.1;
    F = f+0.5*gamma*g.^2;
    surf(X,Y,F)%plot of F
    title('Plot of F=f+0.5*mu*g^2 for mu=0.1')
    subplot(2,2,3)
    gamma = 1;
    F = f+0.5*gamma*g.^2;
    surf(X,Y,F)%plot of F
    title('Plot of F=f+0.5*mu*g^2 for mu=1')
    subplot(2,2,4)
    gamma = 10;
    F = f+0.5*gamma*g.^2;
    surf(X,Y,F)%plot of F
    title('Plot of F=f+0.5*mu*g^2 for mu=10')
    figure(6)
    surf(X,Y,f)
    title('Plot of the objective function f')
    figure(7)
    surf(X,Y,0.5*10*g.^2)
    title('Plot of the penalty function Q_{mu} for mu=10')
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


end