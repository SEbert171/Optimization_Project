clear clc;

import casadi.*
syms X Y
%Settings
%maximum number of NLP iterations:
max_NLP_iterations = 30;
%maximum number of Newton iterations:
max_Newton_iterations = 30;
%test_function_type: ackley and rastrigin are possible
test_function_type='ackley';
%show 3dplot--> very calculation hungry but good to see whats wrong
Show3dplot=true;


% ackley function
if (strcmp(test_function_type, 'ackley'))
    f = -20*exp(-0.2*sqrt(0.5*(X^2+Y^2)))-exp(0.5*(cos(2*pi*X)+cos(2*pi*Y)))+exp(1)+20;
elseif (strcmp(test_function_type, 'rastrigin'))
    f = 20+X^2-10*cos(2*pi*X)+Y^2-10*cos(2*pi*Y);
else
    msg='Test function not recognized. Use ackley or rastrigin.';
    error(msg);
end


% constraint function
g = X^2+Y^2-25;


% starting point
iguess= [-5;-1]; % try -5;-5, -4;-4, -3;-3 ...
x=iguess(1);
y=iguess(2);

% penalty parameter
gamma=10;

% convergence criteria for newton method
e = 10^(-10);

%maximum constraint violation
max_constraint_violation = 10 ^(-4);
tic


% jacobian and hessian derivatives of penalty function F=f+0.5*gamma*(g^2)
[Jp,Hp] = penalty_derivatives(f,g,iguess,gamma);

rguess=iguess; %rguess -> record guesses; iguess -> initial guess

n=0;
k=1;
while (norm(Jp)>10^(-10-k+1) && n<=max_NLP_iterations)
    l=0;
    while (norm(Jp) > e && l<= max_Newton_iterations)
        disp(['Newton iteration: ',num2str(l),' for NLP-iteration ',num2str(n)]);
        [Jp,Hp] = penalty_derivatives(f,g,iguess,gamma);
        %S=inv(Hp); % Search direction--> not needed because one can solve
        %the linear system directly
        iguess = iguess - Hp\Jp;
        rguess=iguess;
        l=l+1;
    end
 constraint_violation = subs(g,[X,Y], [rguess(1),rguess(2)]);
 if (norm(Jp) <= 10^(-10-k+1) && norm(constraint_violation) <= max_constraint_violation)
     break
 end
 n=n+1;
 k=k+1;
 gamma=gamma*10;
 [Jp,Hp] = penalty_derivatives(f,g,rguess,gamma);
end
time_elapsed = toc;
% displaying the results
% rguess(:,n+1)=[];

% plots:
figure(1)
fcontour(f, 'Fill', 'On');
hold on;
plot(rguess(1,:),rguess(2,:),'*-r');
grid on;
fimplicit(g,'r');
colorbar

if(Show3dplot)
figure(2)
fmesh(f)
hold on;
plot3(rguess(1,:),rguess(2,:),subs(f,[X,Y], [rguess(1),rguess(2)]),'*-r');
grid on
fimplicit3([g,0])
colorbar
end

% output:
fprintf('Initial Objective Function Value: %d\n\n',subs(f,[X,Y], [x,y]));
if (norm(Jp) < e)
 fprintf('Minimum succesfully obtained...\n\n');
end

fprintf('Number of Iterations for Convergence: %d\n\n', n);
fprintf('Point of Minima: [%f,%f]\n\n', iguess(1), iguess(2));
fprintf('Objective Function Minimum Value after Optimization: %f\n\n', subs(f,[X,Y], [iguess(1),iguess(2)]));

fprintf('Norm of the constraint violation: %f\n\n', norm(constraint_violation));


disp(['Solution time was ',num2str(time_elapsed),' seconds'])