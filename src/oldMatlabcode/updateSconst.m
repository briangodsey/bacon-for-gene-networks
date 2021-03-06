function[ par ] = updateSconst( D, par, priors, info )

expect = getExpectations( par );

for i = 1:info.nclus

    % calculate the row of the precision for this cluster
    temp1 = 0;
    for tsum = 2:(info.Tmax+1)
        temp1 = temp1 + expect.sig(i,i);
    end
    par.Sconstprec(i) = priors.Sconstprec(i) + temp1;
    %par.Sconstprec
    
    % calculate the mean for this cluster
    temp2 = repmat(0,info.Tmax+1,1);
    for tsum = 2:(info.Tmax+1)
        temp2(tsum) = (par.Fmean(i,tsum)-par.Smean(i,:)*par.Fmean(:,tsum-1)) * expect.sig(i,i);
    end
    sum2 = sum(temp2);
    
    par.Sconstmean(i) =  ( priors.Sconstmean(i)*priors.Sconstprec(i) + sum2 ) / par.Sconstprec(i);
    
end

%par.Sconstmean

end