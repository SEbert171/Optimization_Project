clc;
syms X Y

% ackley function
f = -20*exp(-0.2*sqrt(0.5*(X^2+Y^2)))-exp(0.5*(cos(2*pi*X)+cos(2*pi*Y)))+exp(1)+20;

% constraint function
g = X^2+Y^2-25;

% starting point
iguess= [-1;-1]; % try -5;-5, -4;-4, -3;-3 ...
x=iguess(1);
y=iguess(2);

% penalty parameter
gamma=10^8;

% convergence criteria for newton method
e = 10^(-10);

% sequence of t->0
t = [];
t = [t; 10^(-10)];

i=1;
while t(i)>0
 i=i+1;
 t = [t; 10^(-10-i+1)];
end

% jacobian and hessian derivatives of penalty function F=f+0.5*gamma*(g^2)
[Jp,Hp] = penalty_derivatives(f,g,iguess,gamma);

rguess=iguess; %rguess -> record guesses; iguess -> initial guess

n=0;
k=1;
while (norm(Jp)>t(k))
 while norm(Jp) > e
   [Jp,Hp] = penalty_derivatives(f,g,iguess,gamma);
   S=inv(Hp); % Search direction
   iguess = iguess - S*Jp;
   rguess=iguess;
 end
 
 if norm(Jp) <= t(k)
     break
 end
 
 n=n+1;
 k=k+1;
 gamma=gamma*10;
 [Jp,Hp] = penalty_derivatives(f,g,rguess,gamma);
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