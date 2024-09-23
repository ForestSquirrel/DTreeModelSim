function signal = stimuliChirp(tmax, h, f0, t, f1)
% Function to generate the value of a logarithmic chirp signal at specific timestamps
% signal = stimuliChirp(tmax, h, f0, t, f1)
% tmax - final time (must be positive)
% h - time step (must be positive)
% f0 - starting frequency (must be positive)
% t - timestamp(s) at which to evaluate the signal (must be within [0, tmax])
% f1 - final frequency at tmax (optional)
%     If f1 is not specified, it is set to the Nyquist frequency (1/(2h))

    % Input validation
    if nargin < 4
        error('Not enough input arguments. Usage: stimuliChirp(tmax, h, f0, t, [f1])');
    end

    if tmax <= 0
        error('tmax must be positive.');
    end

    if h <= 0
        error('h must be positive.');
    end

    if f0 <= 0
        error('f0 must be positive.');
    end

    if any(t < 0) || any(t > tmax)
        error('All t values must be within the range [0, tmax].');
    end

    fs = 1 / h;        % Sampling frequency
    f_nyq = fs / 2;    % Nyquist frequency

    % Set f1 to Nyquist frequency if not specified
    if nargin < 5 || isempty(f1)
        f1 = f_nyq;
    else
        if f1 > f_nyq
            error('Specified f1 exceeds Nyquist frequency (%.2f Hz). Decrease f1 or increase h.', f_nyq);
        end
    end

    % Generate the phase of the signal at time t
    if f1 == f0
        % Constant frequency case
        phi = 2 * pi * f0 * t;
    else
        % Logarithmic chirp
        k = log(f1 / f0);
        phi = 2 * pi * f0 * tmax * ((exp(k * t / tmax) - 1) / k);
    end

    % Generate the signal at time t
    signal = sin(phi);
end
