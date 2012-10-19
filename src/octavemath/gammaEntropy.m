function[ entropy ] = gammaEntropy( alpha, beta )

% alpha is the shape
% beta is the inverse scale, like in Gelman
entropy = alpha + log(1/beta) + gammaln(alpha) + (1-alpha)*psi(alpha);

end
