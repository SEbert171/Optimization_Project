function [rguess] = solve_Penalty_NLP_augmented_Newton(F,g, Quadratic_Penalty, gamma,iguess, Newton_terminal_condition, max_Newton_iterations, max_line_search_iterations)
%Solves the Penalty NLP using the augmented Newton step (calculating the
%derivatives using Casadi). F is the penalty function (f+0.5*gamma*g^2) and
%iguess is the initial guess for the Newton iteration. Returns the new
%approximation for a minimum

[Jp,Hp] = calculate_derivatives(F,iguess);
n=1;
while (norm(Jp)> Newton_terminal_condition && n<=max_Newton_iterations)
    disp(['     Newton iteration: ',num2str(n),])
    [Jg,Hg] = calculate_derivatives(g,iguess)
    [JQ, ~] = calculate_derivatives(Quadratic_Penalty,iguess)
    B = [Hp+1/gamma*full(g(iguess(1),iguess(2)))*Hg Jg;Jg.' -gamma];
    c = [-JQ;0];
    solve_linear_system = B\c;
    search_direction = solve_linear_system(1:2);
    t=1;
    k=1;
    while (full(F(iguess(1)+t*search_direction(1),iguess(2)+t*search_direction(2))) >= double(full(F(iguess(1),iguess(2))+0.1*Jp*search_direction))  && k<=max_line_search_iterations )
        t=0.8*t;
        k = k+1;
    end
    disp(['Step length was: ',num2str(t)])
    iguess = iguess + t*search_direction;
    rguess=iguess;
    n=n+1;
end