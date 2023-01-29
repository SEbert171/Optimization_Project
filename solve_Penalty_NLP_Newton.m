function [rguess,constraint_violation] = solve_Penalty_NLP_Newton(f,g,gamma,max_constraint_violation, max_NLP_iterations, iguess, Newton_terminal_condition)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    import casadi.*

    X = casadi.MX.sym('x');
    Y = casadi.MX.sym('y');
    g_eval = Function('g',{[X,Y]},{g});


    rguess=iguess; %rguess -> record guesses; iguess -> initial guess
    [Jp,~] = penalty_derivatives_old(f,g,iguess,gamma);
    n=0;
    k=1;
    %NLP step
    while (n<=max_NLP_iterations)
        l=0;
        %Newton step
        while (norm(Jp) > Newton_terminal_condition && l<= max_Newton_iterations)
            disp(['Newton iteration: ',num2str(l),' for NLP-iteration ',num2str(n)]);
            [Jp,Hp] = penalty_derivatives_old(f,g,iguess,gamma);%still uses the old derivative calculation
            iguess = iguess - Hp\Jp;
            rguess=iguess;
            l=l+1;
        end
        constraint_violation = full(evalf(g_eval(rguess)));
        if (norm(Jp) <= 10^(-10-k+1) && norm(constraint_violation) <= max_constraint_violation)
            break
        end
        n=n+1;
        k=k+1;
        gamma=gamma*10;
    end
end