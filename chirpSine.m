function s= chirpSine(t, f0, f1)
    % Generates a logarithmic chirp sine wave.
    %
    % Usage:
    %   s = chirpSine(t, f0, f1)
    %
    % Inputs:
    %   t  - Time vector (array)
    %   f0 - Initial frequency in Hz
    %   f1 - Final frequency in Hz
    %
    % Output:
    %   s - Logarithmic chirp sine wave
    
    % Ensure t is a column vector
    t = t(:);

    % Initial and final time
    t0 = t(1);
    t1 = t(end);

    % Calculate beta (rate of frequency change)
    beta = log(f1 / f0) / (t1 - t0);

    % Calculate the instantaneous phase
    phi = (2 * pi * f0 / beta) * (exp(beta * (t - t0)) - 1);
    
    % Generate the chirp signal
    s = sin(phi);
end
