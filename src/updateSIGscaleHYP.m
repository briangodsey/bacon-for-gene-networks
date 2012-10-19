function[ par ] = updateSIGscaleHYP( exper, D, par, priors, info )

expect = getExpectations( par.exper{exper} );

par.exper{exper}.SIGscaleHYPshape = priors.exper{exper}.SIGscaleHYPshape ...
    + info.nclus*priors.exper{exper}.SIGshape;
par.exper{exper}.SIGscaleHYPscale = priors.exper{exper}.SIGscaleHYPscale ...
    + sum(sum(expect.SIG));

end