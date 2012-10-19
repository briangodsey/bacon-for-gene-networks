function[ par ] = updateS( D, par, priors, info )

par = updateSmeanMiddleGuy( D, par, priors, info );
par = updateSmeanHYP( D, par, priors, info );
par = updateSmeanMiddleGuy( D, par, priors, info );

par = updateSprecHYP( D, par, priors, info );

expect = getExpectations( par );

for i = 1:info.nclus

    % calculate the row of the precision for this cluster
    temp1 = zeros(info.nclus);
    for tsum = 2:(info.Tmax+1)
        temp1 = temp1 + par.Fmean(:,tsum-1) * expect.SIG(i,i) * par.Fmean(:,tsum-1)';
    end
    par.Sprec{i} = expect.SprecHYP*eye(info.nclus) + temp1;
    
    % calculate the mean for this cluster
    temp2 = repmat(0,info.Tmax+1,info.nclus);
    for tsum = 2:(info.Tmax+1)
        temp2(tsum,:) = (par.Fmean(i,tsum)-par.Sconstmean(i)) * expect.SIG(i,i) * par.Fmean(:,tsum-1)';
    end
    sum2 = sum(temp2,1);

    %par.Smean(i,:) =  ( repmat(par.SmeanHYPmean,1,info.nclus) ...
    %    * (expect.SprecHYP*eye(info.nclus)) ... 
    %    + sum2 ) / par.Sprec{i} ;

    par.Smean(i,:) = ( par.Smean(i,:) * (expect.SprecHYP*eye(info.nclus)) ...
        + sum2 ) / par.Sprec{i} ;
    
end

par = updateSconst( D, par, priors, info );

end