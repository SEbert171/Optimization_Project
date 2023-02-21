clc;

import casadi.*
opti = casadi.Opti();

X = opti.variable();
Y = opti.variable();

gamma=10000;

%f = -20*exp(-0.2*sqrt(0.5*(X^2+Y^2+1e-6)))-exp(0.5*(cos(2*pi*X)+cos(2*pi*Y)))+exp(1)+20; % ackley function
f = 20+X^2-10*cos(2*pi*X)+Y^2-10*cos(2*pi*Y); % rastrigin function
g = X^2+Y^2-25; % constraint function for ackley
% F=f+gamma*g^2;

tic
opti.minimize(f)
 
opti.set_initial(X, -3);
opti.set_initial(Y, 3);

opti.subject_to(g==0)

opti.solver('ipopt');
sol = opti.solve();
 
x=sol.value(X);
y=sol.value(Y);
time_elapsed = toc;

disp('Solution point is:')
disp([x;y]);

disp(['Solution time was ',num2str(time_elapsed),' seconds'])
