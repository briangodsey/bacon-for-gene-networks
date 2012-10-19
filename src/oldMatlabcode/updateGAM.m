function[ par ] = updateGAM( D, par, priors, info )

%expect = getExpectations( par );
expSq = getExpOfSquares( par );

par.gamshape = priors.gamshape + 0.5 * size(info.tuniq,2);

mumat = par.mumean(info.tuniq+1,:);
mumatSq = expSq.mu(info.tuniq+1,:);

Fprecmat = repmat(0, size(par.Fmean) );
for i = 1:size(par.Fmean,2)
    Fprecmat(:,i) = diag(par.Fprec{i});
end
% this next step is complicated; be careful
expSqmshipF = (1./(Fprecmat(:,info.tuniq+1)')*par.mship) + (par.Fmean(:,info.tuniq+1)'*par.mship).^2;

% we have to take only the columns of expect.F which have mu's
%temp = sum(mumatSq  - 2 * mumat .* (par.Fmean(:,info.tuniq+1)'*par.mship) + expSqmshipF);
temp = sum(mumatSq)  ...
    - 2 * sum(mumat .* (par.Fmean(:,info.tuniq+1)'*par.mship)) ...
    + sum(expSqmshipF);

%mumatSq
%mumat
%-2*mumat .* (expect.F(:,info.tuniq+1)'*par.mship)
%expSqmshipF

% abs to correct rounding error for large numbers of clusters
par.gamscale = abs( priors.gamscale + 0.5 * temp );


end