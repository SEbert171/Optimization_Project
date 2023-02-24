function [SOSC_satisfied] = check_SOSC(point,test_function_type, Newton_terminal_condition)
%checks if function_type satisfies the second order sufficient condition
%at point, where the jacobian is small
import casadi.*

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

[~,Hf] = calculate_derivatives(f,point);
eig_Hf = eig(Hf);
negative_eigenvalues = (eig_Hf < Newton_terminal_condition);
if (sum(negative_eigenvalues > 0))
    SOSC_satisfied = false;
else
    SOSC_satisfied = true;
end

end