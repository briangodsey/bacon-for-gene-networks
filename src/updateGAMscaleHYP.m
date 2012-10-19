function[ par ] = updateGAMscaleHYP( exper, D, par, priors, info )

expect = getExpectations( par.exper{exper} );

par.exper{exper}.gamscaleHYPshape = priors.exper{exper}.gamscaleHYPshape ...
    + priors.exper{exper}.gamshape*info.nclus*length(info.Tdat{exper});
par.exper{exper}.gamscaleHYPscale = priors.exper{exper}.gamscaleHYPscale ...
    + sum(sum(expect.gam));

end