function[ par ] = initClust( D, par, info )

tempmumean = repmat(0,size(par.mumean));
for t = info.tuniq
    tnow = info.time == t;
    tind = t+1;
    tempmumean(tind,:) = mean( D(tnow,:), 1 );
end

if( size(tempmumean,2)/info.nclus >10)
    startval = 'cluster';
else
    %startval = 'sample';
    startval = 'uniform';
end

% turn 'onlinephase' to 'off' to speed up
warning off
[idx cent sumd ptd ] = kmeans(tempmumean',info.nclus,'start',startval, ...
    'Replicates',20,'emptyaction','singleton','distance','sqEuclidean', ...
			      'onlinephase','on');
warning on

par.mship(:,:) = 0;
for g = 1:info.G
    par.mship(idx(g),g) = 1;
end

% dims = 2:3;
% clustplot(tempmumean(dims,:)',cent(:,dims),idx);
% dims = 3:4;
% clustplot(tempmumean(dims,:)',cent(:,dims),idx);

%par.mship = mship;

end