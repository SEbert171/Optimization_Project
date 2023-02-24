



theta = linspace(0,2*pi);
x = 5*cos(theta);
y = 5*sin(theta);
if strcmp(test_function_type,'rosenbrock')
    x = sqrt(1.5)*cos(theta);
    y = sqrt(1.5)*sin(theta);
end
if strcmp(test_function_type,'convex')
    x = 1*cos(theta);
    y = 1*sin(theta);
end
hold on
figure(8)
hold on

plot(x,y,'-r', 'LineWidth',1.5)
title('Rosenbrock function and convergence behaviour')
xlabel('X')
ylabel('Y')