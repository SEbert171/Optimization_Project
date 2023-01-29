function [Jf,Hf] = penalty_derivatives_old(f,g,iguess,gamma)
    syms X Y
    % adding the objective and constraint function
    F=f+0.5*gamma*(g^2);

    % calcualting the Jacobian for the function F
    df_dx=diff(F,X); df_dy=diff(F,Y);
    Jf=[df_dx; df_dy];
    Jf=subs(Jf,[X,Y],[iguess(1),iguess(2)]);

    % calculating th Hessian for the function F
    ddf_dxx=diff(df_dx,X); ddf_dxy=diff(df_dx, Y);
    ddf_dyx=diff(df_dy,X); ddf_dyy=diff(df_dy, Y);
    Hf = [ddf_dxx, ddf_dxy; ddf_dyx, ddf_dyy];
    Hf=vpa(subs(Hf, [X,Y],[iguess(1),iguess(2)]),3);
end
