function[ par ] = updateDELTA( t, D, par, priors, info )

% might need to correct for indeterminability
% i.e. set sum of deltas to zero (for each t)

expect = getExpectations( par );

par.deltaprec = priors.deltaprec + 2*info.G*expect.L;

tind = t+1;
if info.N(tind) > 0
    tnow = find(info.time==t);
    trand = tnow(1:end-1);
    tlast = tnow(end);

    for i = trand
        sum1 = sum( D(i,:) - par.mumean(tind,:) );
        sum2 =  sum( -D(tlast,:)) + sum(par.mumean(tind,:)) - info.G*sum(par.deltamean( setdiff(trand,i) ));
        par.deltamean(i) = (priors.deltamean(i) * priors.deltaprec(i) + expect.L*(sum1+sum2)) ./ par.deltaprec(i);
    end
    
    par.deltamean(tlast) = -sum(par.deltamean(trand));
end

end