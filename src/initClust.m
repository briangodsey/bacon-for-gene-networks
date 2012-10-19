function[ par ] = initClust( D, par, info )

tempmumean = repmat(0,size(par.mumean));
tempvar = repmat(0,size(par.mumean));
for t = info.tuniq
    tnow = info.time == t;
    tind = t+1;
    tempmumean(tind,:) = mean( D(tnow,:), 1 );
    tempvar(tind,:) = var( D(tnow,:), 1 );
end
%%1./tempvar

if( size(tempmumean,2)/info.nclus >10)
    startval = 'cluster';
else
    %startval = 'sample';
    startval = 'uniform';
end

%% % turn 'onlinephase' to 'off' to speed up
%% warning off
%% [idx cent sumd ptd ] = kmeans(tempmumean',info.nclus,'start',startval, ...
%%     'Replicates',20,'emptyaction','singleton','distance','sqEuclidean', ...
%% 			      'onlinephase','on');
%% warning on

%par.mship(:,:) = 0;
%for g = 1:info.G
%  %%par.mship(idx(g),g) = 1;
%  %%mod(g,size(par.mship,1)+1)
%  par.mship(mod(g,size(par.mship,1))+1,g) = 1;
%end

par.mship = rand(info.nclus,info.G) > (1-1/info.nclus);
mship = par.mship*1;
par.mship(:,:) = 0;

for i = 1:30
%%while sum(sum((par.mship-mship).^2)) > 1e-10 
  par.mship = mship;

  majormembs = sum(mship>0.5,2);
  minormembs = sum(mship>0.1,2);
  if any(majormembs==0)
    emptyclusts = find(majormembs==0);
    bigclusts = find(minormembs==max(minormembs));
    ugenes = find(mship(bigclusts(1),:));
    mship(emptyclusts(1),ugenes(1:floor(length(ugenes)/2))) = 1;
    mship(bigclusts(1),ugenes(1:floor(length(ugenes)/2))) = 0;
  end

  for g = 1:info.G
    dpt = tempmumean(:,g);
    mns = tempmumean*(mship./repmat(sum(mship,2),1,info.G))';
    dists = sum( (repmat(dpt,1,info.nclus) - mns).^2, 1 );
    meandist = mean(std(D))/10; %mean(dists)^0.5;
    vec = meandist ./ max(dists,1e-8);
    mship(:,g) = vec/sum(vec);
    %% lowestdist = find(dists==min(dists));
    %% mship(lowestdist(1),g) = 1;
  end

  %%mship
end

par.mship = mship;

% dims = 2:3;
% clustplot(tempmumean(dims,:)',cent(:,dims),idx);
% dims = 3:4;
% clustplot(tempmumean(dims,:)',cent(:,dims),idx);

%par.mship = mship;

end