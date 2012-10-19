function[ par llh iter ] = vbmodelfit( D, priors, info )
%function[ par ] = vbmodelfit( D, priors, info )

% this parameter sets the number of initial updates to be done, as well as
% the number of iterations between llh calculation
initialIterations = 1;
displayIteration = 10;

par = priors;

if( info.G ~= size(D,2) )
    disp('Adjusting info.G to reflect data set size.\n');
    info.G = size(D,2);
end

par = initClust( D, par, info );

% initialize some parameters to avoid possible problems
for i = 1:initialIterations
    for t = [0 info.tuniq]
        par = updateMU( t, D, par, priors, info );
        par = updateDELTA( t, D, par, priors, info );
        par = updateGAM( D, par, priors, info );
        par = updateF( t, D, par, priors, info );
    end

    par = updateL( D, par, priors, info );
    par = updateS( D, par, priors, info );
    par = updateSIG( D, par, priors, info );
    %par = getMembership( par, info );

    for t = [fliplr(info.tuniq) 0]
        par = updateMU( t, D, par, priors, info );
        par = updateDELTA( t, D, par, priors, info );
        par = updateGAM( D, par, priors, info );
        par = updateF( t, D, par, priors, info );
    end
    par = updateL( D, par, priors, info );
    par = updateS( D, par, priors, info );
    par = updateSIG( D, par, priors, info );
end

par = getMembership( par, info );

% for i=1:10
%     par = updateVB( D, par, priors, info );
%     llh = geneInterLikelihood( D, par, priors, info );
% end

% initialize log likelihood
llh = geneInterLikelihood( D, par, priors, info );

iter = 0;
oldllh = llh(1)-10^15;
%for i = 1:1000
while( abs((llh-oldllh)/llh) > 10^-8 )
    par = updateVB( D, par, priors, info );
    iter = iter+1;

    oldllh = llh;
    llh = geneInterLikelihood( D, par, priors, info );
    llhdiff = abs((llh-oldllh)/llh);

    if( llhdiff < -10^3 )
        llh-oldllh
        iter
        break;
    end
    
    if( iter > 1000 )
        break;
    end
        
    if round(iter/displayIteration)==iter/displayIteration
        iter
    end
end

totaliter = iter

end
