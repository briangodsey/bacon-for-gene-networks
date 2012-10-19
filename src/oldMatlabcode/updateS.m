function[ par ] = updateS( D, par, priors, info )

% SEE FUNCTION DEFINITION FOR updateSprecalc FOR MORE INFO
% not a good idea to use this because of singularities
%par = updateSprecalc( D, par, priors, info );

expect = getExpectations( par );

for i = 1:info.nclus

    % calculate the row of the precision for this cluster
    temp1 = repmat(0,info.nclus,info.nclus);
    for tsum = 2:(info.Tmax+1)
        temp1 = temp1 + par.Fmean(:,tsum-1) * expect.sig(i,i) * par.Fmean(:,tsum-1)';
    end
    par.Sprec{i} = priors.Sprec{i} + temp1;
    %par.Sprec{i}
    
    % calculate the mean for this cluster
    temp2 = repmat(0,info.Tmax+1,info.nclus);
    for tsum = 2:(info.Tmax+1)
        %par.Fmean(i,tsum)
        %par.Sconstmean(i)
        %expect.sig(i,i)
        %par.Fmean(:,tsum-1)
        temp2(tsum,:) = (par.Fmean(i,tsum)-par.Sconstmean(i)) * expect.sig(i,i) * par.Fmean(:,tsum-1)';
    end
    sum2 = sum(temp2,1);
    
    %i
    %temp2
    %sum2
    
    %priors.Smean(i,:)*priors.Sprec{i} + sum2;
    
    par.Smean(i,:) =  ( priors.Smean(i,:)*priors.Sprec{i} + sum2 ) / par.Sprec{i} ;
    
end

par = updateSconst( D, par, priors, info );

%par.Smean
%par.Sconstmean


end