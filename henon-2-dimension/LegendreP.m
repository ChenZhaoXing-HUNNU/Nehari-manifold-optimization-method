function P = legendreP(n, x)
    % Compute the value n-Legendre in x
    P = legendre(n, x);
    P = P(1, :);
end