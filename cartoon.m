% Filter parameters
N = 1536; % Number of points in PFB window
K = 128;  % Number of points in resulting spectrum
ntaps = 4;
tapsize = N/ntaps;

% The sinc function
function f = filter(N, K)
    n = [0:N-1]';

    % The sinc function
    x = pi*(n+1-N/2)/K;
    s = sin(x)./x;
    s(N/2) = 1;

    % The window function (Hanning)
    h = sin(pi*(n+1)/(N+1)).^2;

    % The final filter
    f = s.*h;
end

% The signal function
function x = signal(N)
    n = [0:N-1]';
    x = cos(0.01*n) + randn(N,1)/2.0;
end

% Do a PFB!
sig = signal(N);
fil = filter(N, K);
filsig = fil .* sig;

summed = sum(reshape(filsig, [tapsize, ntaps]), 2);
ffted = abs(fft(summed));

data = [sig, fil, filsig];

result = [summed, ffted];

save "cartoon.dat" data result
