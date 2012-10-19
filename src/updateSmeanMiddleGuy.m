function[ par ] = updateSmeanMiddleGuy( D, par, priors, info )

expect = getExpectations( par );

SprecMAT = zeros(info.nclus,info.nclus);
for i = 1:info.nclus
    SprecMAT(i,:) = diag(par.Sprec{i})';
end
par.SmeanMiddleGuyprec = expect.SmeanHYPdoublePREC + SprecMAT;

par.SmeanMiddleGuymean = ( expect.SmeanHYPdoublePREC*par.SmeanHYPmean ...
    + SprecMAT.*par.Smean ) ./ par.SmeanMiddleGuyprec;
end