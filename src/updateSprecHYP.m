function[ par ] = updateSprecHYP( D, par, priors, info )

%expect = getExpectations( par );

par.SprecHYPshape = priors.SprecHYPshape + 0.5*info.nclus^2;

% SprecMat = zeros(info.nclus);
% for i = 1:info.nclus
%     SprecMat(i,:) = diag(par.Sprec{i});
% end
expSsquared = sum(sum(par.SmeanMiddleGuymean.^2 + 1./par.SmeanMiddleGuyprec));
expSxShyp = sum(sum(par.SmeanMiddleGuymean * par.SmeanHYPmean));
expSmeanHYPsq = sum(sum(par.SmeanHYPmean^2 + 1/par.SmeanHYPprec));

par.SprecHYPscale = priors.SprecHYPscale + 0.5 ...
    * ( expSsquared -2*expSxShyp + expSmeanHYPsq );

end
