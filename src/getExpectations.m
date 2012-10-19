function[ expect ] = getExpectations( par )

%% expect.SprecHYP = par.SprecHYPshape / par.SprecHYPscale;
%% 
%% expect.SmeanHYPdoublePREC = par.SmeanHYPdoublePRECshape / par.SmeanHYPdoublePRECscale;

expect.gamscaleHYP = par.gamscaleHYPshape ./ par.gamscaleHYPscale;
expect.gam = par.gamshape ./ par.gamscale;

expect.L = par.Lshape/par.Lscale;

expect.SIGscaleHYP = par.SIGscaleHYPshape ./ par.SIGscaleHYPscale;
expect.SIG = diag(par.SIGshape ./ par.SIGscale);
%expect.SIG(expect.SIG<1e-10 & expect.SIG>0) = 1e-2;
%expect.SIG = eye(size(par.SIGshape,2));

end