clc;

import casadi.*

opti = casadi.Opti();

X = opti.variable();
Y = opti.variable();

% ackley function
f = -20*exp(-0.2*sqrt(0.5*(X^2+Y^2)))-exp(0.5*(cos(2*pi*X)+cos(2*pi*Y)))+exp(1)+20;

% constraint function
g = X^2+Y^2-25;

% starting point
iguess= [-1;-1]; % try -5;-5, -4;-4, -3;-3 ...
x=iguess(1);
y=iguess(2);

% penalty parameter
gamma=10;

% sequence of t->0
t = [];
t = [t; 10^(-10)];

i=1;
while t(i)>0
 i=i+1;
 t = [t; 10^(-10-i+1)];
end

% jacobian and hessian derivatives of penalty function F=f+0.5*gamma*(g^2)
[Jp] = penalty_derivatives(f,g,iguess,gamma);

n=0;
k=1;
while (norm(Jp)>t(k))

 opti.minimize(f+0.5*gamma*(g^2))
 opti.solver('ipopt');
 sol = opti.solve();
 
 x=sol.value(X);
 y=sol.value(Y);
 iguess=[x;y];
 
 [Jp] = penalty_derivatives(f,g,iguess,gamma);
 
 if (norm(Jp) <= t(k) || n==15)
     break
 end
 
 n=n+1;
 k=k+1;
 gamma=gamma*10;
 [Jp] = penalty_derivatives(f,g,iguess,gamma);
 
 opti.minimize(f+0.5*gamma*(g^2))
 opti.solver('ipopt');
 sol = opti.solve();
 
 x=sol.value(X);
 y=sol.value(Y);
 iguess=[x;y];
 
 [Jp] = penalty_derivatives(f,g,iguess,gamma);
end

% displaying the results
% rguess(:,n+1)=[];

% plots:
fcontour(f, 'Fill', 'On');
hold on;
plot(rguess(1,:),rguess(2,:),'*-r');
grid on;

% output:
fprintf('Initial Objective Function Value: %d\n\n',subs(f,[X,Y], [x,y]));
if (norm(Jp) < e)
 fprintf('Minimum succesfully obtained...\n\n');
end

fprintf('Number of Iterations for Convergence: %d\n\n', n);
fprintf('Point of Minima: [%f,%f]\n\n', iguess(1), iguess(2));
fprintf('Objective Function Minimum Value after Optimization: %f\n\n', subs(f,[X,Y], [iguess(1),iguess(2)]));