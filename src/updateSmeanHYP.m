function[ par ] = updateSmeanHYP( D, par, priors, info )

expect = getExpectations( par );

par.SmeanHYPprec = expect.SmeanHYPdoublePREC * (1 + info.nclus^2);
par.SmeanHYPmean = (expect.SmeanHYPdoublePREC*priors.SmeanHYPmean ...
    + sum(sum(expect.SmeanHYPdoublePREC * par.SmeanMiddleGuymean)) ) ...
    / par.SmeanHYPprec;

end