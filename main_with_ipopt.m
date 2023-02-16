clc;

import casadi.*
opti = casadi.Opti();

X = opti.variable();
Y = opti.variable();
z = [X;Y]; 

%Settings
%maximum number of NLP iterations:
max_NLP_iterations = 30;
%test_function_type: ackley and rastrigin are possible
test_function_type='ackley';
%show 3dplot--> very calculation hungry but good to see whats wrong
Show3dplot=true;

if (strcmp(test_function_type, 'ackley'))
    f = -20*exp(-0.2*sqrt(0.5*(X^2+Y^2 +10^(-3) )))-exp(0.5*(cos(2*pi*X)+cos(2*pi*Y)))+exp(1)+20; % ackley function
    g = X^2+Y^2-25; % constraint function for ackley
elseif (strcmp(test_function_type, 'rastrigin'))
    f = 20+X^2-10*cos(2*pi*X)+Y^2-10*cos(2*pi*Y); % rastrigin function
    g = X^2+Y^2-26.2144; % constraint function for rastrigin
else
    msg='Test function not recognized. Use ackley or rastrigin.';
    error(msg);
end

% starting point
iguess= [4;4]; % try -5;-5, -4;-4, -3;-3 ...
x=iguess(1);
y=iguess(2);

% penalty parameter
gamma0=1;

tic
% jacobian derivative of penalty function F=f+0.5*gamma*(g^2)
[Jp] = penalty_derivatives(iguess,gamma0);

k=1;
while (norm(Jp)>10^(-10-k+1) && k<=max_NLP_iterations)
 
 opti.minimize(f+0.5*gamma0*(g^2))
 opti.solver('ipopt');
 sol = opti.solve();
 
 x=sol.value(X);
 y=sol.value(Y);
 iguess=[x;y];
 gamma0=gamma0*10;
 
 [Jp] = penalty_derivatives(iguess,gamma0);

 if (norm(Jp) <= 10^(-10-k+1) || k==15)
     break
 end
 k=k+1;
end

time_elapsed = toc;


% % displaying the results
% 
% % plots:
% figure(1)
% fcontour(f, 'Fill', 'On');
% hold on;
% plot(rguess(1,:),rguess(2,:),'*-r');
% grid on;
% fimplicit(g,'r');
% colorbar
% 
% if(Show3dplot)
% figure(2)
% fmesh(f)
% hold on;
% plot3(iguess(1,:),iguess(2,:),subs(f,[X,Y], [iguess(1),iguess(2)]),'*-r');
% grid on
% fimplicit3([g,0])
% colorbar
% end
% 
% output:
% fprintf('Initial Objective Function Value: %d\n\n',subs(f,[X,Y], [x,y]));

fprintf('Number of Iterations for Convergence: %d\n\n', k);
fprintf('Point of Minima: [%f,%f]\n\n', iguess(1), iguess(2));
% fprintf('Objective Function Minimum Value after Optimization: %f\n\n', subs(f,[X,Y], [iguess(1),iguess(2)]));

disp(['Solution time was ',num2str(time_elapsed),' seconds']);