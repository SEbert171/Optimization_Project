function [Jf,Hf] = calculate_derivatives(f, x0)
%calculates the jacobian for a 2-dimensional Casadi function f at the
%point x0 in R^2
import casadi.*
 




 %jacobian of the function f
 JF = jacobian(f);
 %Hessian of the function f, i.e jacobian of the jacobian of f
 HF = jacobian(JF);
 % calcualting the Jacobian for the function f
 Jf = full(JF(x0(1),x0(2),2));
 Hessian = full(HF(x0(1),x0(2),1,1));
 Hf = Hessian(:,1:2);
end