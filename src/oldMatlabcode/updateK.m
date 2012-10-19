function[ par ] = updateK( D, par, priors, info )

expect = getExpectations( par );

par.Kprec = priors.Kprec + expect.gam .* sum(par.Fmean'*par.mship,1).^2;

temp = sum( expect.mu .* (expect.F'*par.mship) ,1 );
par.Kmean = (priors.Kprec.*priors.Kmean + expect.gam .* temp) ./ par.Kprec;

end