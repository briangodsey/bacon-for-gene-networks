function[ par ] = stepClustering( par, info )

par = getMembership( par, info );

par.Fmean = par.mship' \ par.mumean';
%par.Fmean(1,:) = 1;

end
