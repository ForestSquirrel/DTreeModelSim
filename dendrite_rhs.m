function [dUdt, dVdt] = dendrite_rhs(u, v, params, coupling_term, Stim)
    % RHS for dendrite
    alpha = params.alpha;
    b = params.b;
    Tau = params.Tau;
    TauR = params.TauR;
    NaX = params.NaX;
    dUdt = Tau * (NaX *(u*(u - 1)*(1 - alpha*u) - v) + coupling_term + Stim);
    dVdt = TauR * b*u;
end
