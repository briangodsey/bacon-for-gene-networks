function[ par ] = updateVB( D, par, priors, info )

%llh = [];
%llh(1) = geneInterLikelihood( D, par, priors, info );

par = updateL( D, par, priors, info );
par = updateGAM( D, par, priors, info );
par = updateSIG( D, par, priors, info );
par = updateS( D, par, priors, info );
%par = updateK( D, par, priors, info );
par = getMembership( par, info );

par = forwardPass( D, par, priors, info );

par = updateL( D, par, priors, info );
par = updateGAM( D, par, priors, info );
par = updateSIG( D, par, priors, info );
par = updateS( D, par, priors, info );
%par = updateK( D, par, priors, info );
par = getMembership( par, info );

par = backwardPass( D, par, priors, info );

par = updateL( D, par, priors, info );
par = updateGAM( D, par, priors, info );
par = updateSIG( D, par, priors, info );
par = updateS( D, par, priors, info );
%par = updateK( D, par, priors, info );


% if( any(llh(2:17)-llh(1:16)<-10^-2) )
%     llh(2:17)-llh(1:16)
%     llh
% end


end