function[ expect ] = getExpectations( par )

expect.gam = par.gamshape ./ par.gamscale;
expect.L = par.Lshape/par.Lscale;
expect.sig = diag(par.SIGshape ./ par.SIGscale);

end