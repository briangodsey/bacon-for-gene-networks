function[ par ] = updateSIG( D, par, priors, info )

par.SIGshape = priors.SIGshape + repmat(size(info.tuniq,2),1,info.nclus)/2;

tempsum = repmat(0,1,info.nclus);
for i = 1:info.nclus
    temp = repmat(0,6,info.Tmax);
    for t = 1:info.Tmax
        tind = t+1;
        expFsq = par.Fmean(i,tind)^2 + 1/par.Fprec{tind}(i,i);
        
        %expSFsq1 = ( par.Smean(i,:)*par.Fmean(:,tind-1) )^2 + 1./(abs(par.Smean(i,:))*diag(par.Fprec{tind-1}));
        %expSFsq2 = (par.Smean(i,:)*par.Fmean(:,tind-1))^2 + ( ( par.Smean(i,:).^2+1./diag(par.Sprec{i})' ) * ( par.Fmean(:,tind-1).^2+1./(diag(par.Fprec{tind-1})) ) - par.Smean(i,:).^2*par.Fmean(:,tind-1).^2 );

        % this calculates the expectation of (s1*f1+s2*f2+...)^2
        expSsq = par.Smean(i,:).^2 + 1./diag(par.Sprec{i})';
        expFminus1sq = par.Fmean(:,tind-1).^2 + 1./diag(par.Fprec{tind-1});
        SFelementwise = par.Smean(i,:).*par.Fmean(:,tind-1)';
        expSFsq = expSsq*expFminus1sq + sum(sum(SFelementwise'*SFelementwise)) - trace(SFelementwise'*SFelementwise);

        %expSsq*expFminus1sq
        %sum(sum(SFelementwise'*SFelementwise))
%         sum(- trace(SFelementwise'*SFelementwise))
%         SFelementwise*SFelementwise'
%         sum((par.Fmean(:,tind-1).^2)' .* par.Smean(i,:).^2)
%         expSsq*expFminus1sq
%         par.Smean
%         par.Smean(i,:) 
%         par.Smean(i,:).^2 
%         1./diag(par.Sprec{i})'
%         par.Fmean(:,tind-1).^2 
%         1./diag(par.Fprec{tind-1})
%         
        expSconstSq = par.Sconstmean(i)^2 + 1./par.Sconstprec(i);
        expFtSF = par.Fmean(i,tind) * par.Smean(i,:)*par.Fmean(:,tind-1);
        expFtSconst = par.Fmean(i,tind) * par.Sconstmean(i);
        expSFSconst = par.Smean(i,:)*par.Fmean(:,tind-1) * par.Sconstmean(i);
        %temp(t) = expFsq +expSFsq +expSconstSq -2*expFtSF -2*expFtSconst +2*expSFSconst;
        
        % save the elements of the sum into a column for each t; adding up
        % similar numbers first helps avoid rounding errors; we later sum
        % the columns first and then the rows
        temp(:,t) = [ expFsq; expSFsq; expSconstSq; -2*expFtSF; -2*expFtSconst; 2*expSFSconst ];
        %sum(temp,1)
    end
        
    tempsum(i) = sum(sum(temp,2));
end
par.SIGscale = priors.SIGscale + 0.5 * tempsum;

end