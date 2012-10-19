function[ par ] = updateL( exper, D, par, priors, info )

%expect = getExpectations( par );
expSq = getExpOfSquares( par.exper{exper} );

par.exper{exper}.Lshape = priors.exper{exper}.Lshape + 0.5 * info.G * sum(sum(info.N{exper}));

mumat = par.exper{exper}.mumean(info.timeD{exper},:);
mumatSq = expSq.mu(info.timeD{exper},:);

temp = sum(sum(D{exper}.^2)) + sum(sum(mumatSq)) ...
    - 2 * sum(sum(D{exper}.*mumat));

% abs to correct rounding error for large numbers of clusters
par.exper{exper}.Lscale = abs( priors.exper{exper}.Lscale + 0.5 * temp );

end