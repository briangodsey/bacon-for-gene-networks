function[ par ] = backwardPass( D, par, priors, info )

for t = fliplr( [0:info.Tmax] )
    if( sum(info.time == t) >= 1 )
        par =  updateDELTA( t, D, par, priors, info );
        par =  updateMU( t, D, par, priors, info );
    end
    par =  updateF( t, D, par, priors, info );
end

end