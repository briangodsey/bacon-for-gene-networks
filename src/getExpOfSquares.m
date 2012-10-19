function[ expSq ] = getExpOfSquares( par )

expSq.mu = 1./par.muprec + par.mumean.^2;

end