function[ expSq ] = getExpOfSquares( par )

expSq.mu = 1./par.muprec + par.mumean.^2;
expSq.delta = 1./par.deltaprec + par.deltamean.^2;

end