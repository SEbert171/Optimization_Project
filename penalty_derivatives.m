function [Jf] = penalty_derivatives(iguess,gamma0)
 import casadi.*
 X = MX.sym('X',1);
 Y = MX.sym('Y',1);
 
 % optimize w.r.t. z
 z = [X;Y]; 
 
 % parameter – penalty parameter
 gamma = MX.sym('gamma',1);
 
 % ackley function
 f = -20*exp(-0.2*sqrt(0.5*(X^2+Y^2 +10^(-3)   )))-exp(0.5*(cos(2*pi*X)+cos(2*pi*Y)))+exp(1)+20;
 
 % constraint function for ackley
 g = X^2+Y^2-25;
 
 % adding the objective and constraint function
 F=f+0.5*gamma*(g^2);
 JF = F.jacobian(z);
 
 % create casadi function
 JF_fun = Function('JF_fun',{z,gamma},{JF});

 % calcualting the Jacobian for the function F
 Jf=full(JF_fun(iguess,gamma0));
end