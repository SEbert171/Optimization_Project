function [guess, constraint_violation] = Ipopt_solve_NLPs(f,g,gamma,max_constraint_violation, max_NLP_iterations, iguess)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    import casadi.*
    %opti = casadi.Opti();

    X = casadi.MX.sym('x');
    Y = casadi.MX.sym('y');
    
    g_eval = Function('g',{[X,Y]},{g});

    n=0;
    k=1;
    while true

        Penalty_NLP = struct('x',[X,Y],'f',f+0.5*gamma*(g^2));
        S = nlpsol('S','ipopt',Penalty_NLP);
        r = S('x0',iguess);
        guess = full(evalf(r.x));
        %opti.minimize(f+0.5*gamma*(g^2))
        %opti.solver('ipopt');
        %sol = opti.solve();
 
        %disp(sol)
        %x=sol.value(X);
        %y=sol.value(Y);
        %guess=[x;y];
        constraint_violation = full(evalf(g_eval(guess)));
        if (norm(constraint_violation) <= max_constraint_violation || n == max_NLP_iterations)
            break
        end
        n=n+1;
        k=k+1;
        gamma=gamma*100;
    end
end