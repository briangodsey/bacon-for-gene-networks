function[ par ] = updateSprecalc( D, par, priors, info )
% WARNING: THIS FUNCTION HAD SOME PROBLEMS WITH NEAR-SINGULAR MATRICES, BUT
% EITHER WAY THIS FUNCTION WAS INTENDED ONLY TO SPEED UP CALCULATIONS BY
% SOLVING FOR THE CONSTANT (Sconst) AND INTERACTION (S) TERMS
% SIMULTANEOUSLY. USE AT YOUR OWN RISK

expect = getExpectations( par );

for i = 1:info.nclus

    % calculate the row of the precision for this cluster
    temp1 = repmat(0,info.nclus+1,info.nclus+1);
    for tsum = 2:(info.Tmax+1)
        FmeanAug = [ 1; par.Fmean(:,tsum-1) ];
        temp1 = temp1 + FmeanAug * expect.sig(i,i) * FmeanAug';
    end
    sprec = temp1;
    
    % calculate the mean for this cluster
    temp2 = repmat(0,info.Tmax+1,info.nclus+1);
    for tsum = 2:(info.Tmax+1)
        FmeanAugBefore = [ 1; par.Fmean(:,tsum-1) ];
        FmeanAfter =  par.Fmean(i,tsum);
        temp2(tsum,:) = FmeanAfter * expect.sig(i,i) * FmeanAugBefore';
    end
    sum2 = sum(temp2,1);
    
    %i
    %temp1
    %temp2

    %priors.Smean(i+1,:)*priors.Sprec{i+1} + sum2
    %inv(par.Sprec{i+1})
    %par.Sprec{i+1}
    
    smean = sum2 / sprec;
    
    par.Smean(i,:) = smean(2:end);
    par.Sconstmean(i) = smean(1);
    
end

end