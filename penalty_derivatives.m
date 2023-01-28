function [Jf] = penalty_derivatives(f,g,iguess,gamma)

 import casadi.*

 opti = casadi.Opti();

 X = opti.variable();
 Y = opti.variable();

 % adding the objective and constraint function
 F=f+0.5*gamma*(g^2);
    
 % calcualting the Jacobian for the function F at iguess
 Jf=jacobian(F(X,Y),iguess);
end