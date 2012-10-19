function[ par ] = updateL( D, par, priors, info )

%expect = getExpectations( par );
expSq = getExpOfSquares( par );

par.Lshape = priors.Lshape + 0.5 * info.G * sum(sum(info.N));

mumat = par.mumean(info.time+1,:);
deltamat = repmat( par.deltamean', 1, info.G );
mumatSq = expSq.mu(info.time+1,:);
deltamatSq = repmat( expSq.delta', 1, info.G );

%temp1 = sum(sum(D.^2 + deltamatSq + mumatSq - 2 * D.*deltamat - 2 * D.*mumat + 2 * deltamat.*mumat))

temp = sum(sum(D.^2)) + sum(sum(deltamatSq)) ...
    + sum(sum(mumatSq)) - 2 * sum(sum(D.*deltamat)) ...
    - 2 * sum(sum(D.*mumat)) + 2 * sum(sum(deltamat.*mumat));

% abs to correct rounding error for large numbers of clusters
par.Lscale = abs( priors.Lscale + 0.5 * temp );

end