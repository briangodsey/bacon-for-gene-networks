function[ par ] = getMembership( par, info )

mship = zeros( info.nclus, info.G );

%% soft Bayesian clustering
logtemp = zeros(info.nclus,info.G);
%logAmat = [];
%logDmat = [];
%distmat = [];
for g = 1:info.G
    dpt = par.mumean(info.tuniq+1,g);
    dptprec = diag( par.muprec(info.tuniq+1,g) );

    for i = 1:info.nclus
        ctr = par.Fmean(i,info.tuniq+1)';

        %% use the precision of the genes, par.muprec
        %A = (2*pi)^(-0.5*size(info.tuniq,2)) * det(dptprec)^0.5;
        %D = exp( -0.5*(dpt - ctr)' * dptprec * (dpt - ctr) );

        logA = (-0.5*size(info.tuniq,2))*log(2*pi) + 0.5*logdet(dptprec);
        logD = -0.5 * (dpt - ctr)' * dptprec * (dpt - ctr);
        
        %distmat(i,g) = sum((dpt-ctr).^2)
        %logAmat(i,g) = logA
        %logDmat(i,g) = logD
        logtemp(i,g) = logA + logD;
        %temp(i,g) = exp(logtemp(i,g))
    end
end
%tempsum = sum(temp,1);
%tempsummat = repmat(tempsum,info.nclus,1);
%mship = temp./tempsummat;

for i = 1:info.nclus
    logrow = repmat(logtemp(i,:),info.nclus,1);
    mship(i,:) = 1 ./ sum(exp(logtemp-logrow),1);
end

emptyclusters = find( sum(mship,2)==0 )';
if( size(emptyclusters,1) ~= 0 )
    maxlogtemp = max(logtemp,[],1);
    for i = emptyclusters
        loneliestgene = find( maxlogtemp == min(maxlogtemp) );
        mship(:,loneliestgene) = 0;
        mship(i,loneliestgene) = 1;
        maxlogtemp(loneliestgene) = max(maxlogtemp)+1;
    end
end
maxmembership = max(sum(mship,1));

% % correct for the problem where all likelihoods are too small so that they
% % sum to zero, giving improper membership values
% for i = find(tempsum==0)
%     % assign equal membership to all clusters
%     mship(:,i) = 1/info.nclus;
% end

%% hard (Bayesian) clustering (based on likelihood, not distance)
% for i = 1:info.G
%     % or assign to the cluster with the highest likelihood
%     % needs logA, logD, and logtemp from above
%     mship(:,i) = logtemp(:,i)==max(logtemp(:,i));
% end

%% hard clustering with kmeans (not working well)

% for g = 1:info.G
%     %for j = 1:info.nclus
%     [ M I ] = min( sum( (repmat(par.mumean(:,g),1,info.nclus) - par.Fmean(:,info.tuniq+1)').^2 ) );
%     %end
%     mship(I,g) = 1;
% end

% [idx cent sumd ptd ] = kmeans(par.mumean(info.tuniq,:)',info.nclus,'start',par.Fmean(:,info.tuniq+1),'emptyaction','singleton','distance','sqEuclidean');
% for g = 1:info.G
%     mship(idx(g),g) = 1;
% end

%% end
par.mship = mship;

end