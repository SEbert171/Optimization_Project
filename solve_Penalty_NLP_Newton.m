function [rguess] = solve_Penalty_NLP_Newton(F, Q, iguess, Newton_terminal_condition, max_Newton_iterations, max_line_search_iterations)
%Solves the Penalty NLP using the exact Newton step (calculating the
%derivatives using Casadi). F is the penalty function (f+0.5*gamma*g^2) and
%iguess is the initial guess for the Newton iteration. Returns the new
%approximation for a minimum

starting_point = iguess;
[Jp,Hp] = calculate_derivatives(F,iguess);
n=1;
while (norm(Jp)> Newton_terminal_condition && n<=max_Newton_iterations)
    %disp(['    Newton iteration: ',num2str(n)]);
    search_direction = (-Hp\Jp.');
    %Line search with Armijo Backtracking
    t=1;
    k=1;
    while (full(F(iguess(1)+t*search_direction(1),iguess(2)+t*search_direction(2))) >= double(full(F(iguess(1),iguess(2))+0.1*Jp*search_direction))  && k<=max_line_search_iterations )
          t=0.8*t;
          k = k+1;
    end
    %disp(['         Step length was: ',num2str(t)])
    iguess = iguess + t.*search_direction;
    if (full(Q(iguess(1),iguess(2))) > full(Q(starting_point(1),starting_point(2))))
        %disp('Newton step converged away from feasible set! Return.')
        rguess = starting_point;
        break
    end   
    rguess=iguess;
    [Jp,Hp] = calculate_derivatives(F,iguess);
    
    n=n+1;
end
