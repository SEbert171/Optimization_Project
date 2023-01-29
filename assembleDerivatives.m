function [Jf,Hf] = assembleDerivatives (f,g,iguess,gamma)
%UNTITLED2 Assembles derivatives with the help of Casadis Jacobian and
%Hessian Function
import casadi.*

 

 X = casadi.MX.sym('x');
 Y = casadi.MX.sym('y');

 % adding the objective and constraint function
 F=f+0.5*gamma*(g^2);
    
 % calcualting the Jacobian for the function F at iguess
 Jf_sym=Function('jacobian',{[X,Y]},{jacobian(F,[X,Y])});
 Hf_sym = Function('hessian',{[X,Y]},{hessian(F,[X,Y])});
 Jf = full(Jf_sym(iguess));
 Hf = full(Hf_sym(iguess));
end