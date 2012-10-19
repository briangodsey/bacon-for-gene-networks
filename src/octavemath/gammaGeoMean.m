function[ geomean ] = gammaGeoMean( alpha, beta )

% alpha is the shape
% beta is the inverse scale, like in Gelman
geomean = psi(alpha) - log(beta);

end
