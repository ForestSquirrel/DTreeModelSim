function [dUdt, dVdt] = soma_rhs(u, v, params, coupling_term)
    % RHS for soma
    alpha = params.alpha;
    b = params.b;
    Tau = params.Tau;
    TauR = params.TauR;
    dUdt = Tau * (u*(u - 1)*(1 - alpha*u) - v + coupling_term);
    dVdt = TauR * b*u;
end