function[ auroc aupr inferredints ] = checkInteractionResults( par, goldstd )

  %% evaluate sensitivity and specificity
  simints = goldstd; %%getInteractionSignificance( goldstd, 0, eye(size(par.exper{1}.mumean,2)) );
  [ inferredints Ssigs ] = getInteractionSignificance( par.Smean, par.Sprec, par.mship );

  %% remove self-interactions
  for clus = 1:size(inferredints,1)
    inferredints(clus,clus) = 0;
  end

  infintvec = abs( reshape(inferredints,[],1) );
  [ sortints intorder ] = sort(infintvec,"descend");

  simintvec = reshape(simints,[],1);
  sortsimints = simintvec(intorder,:);

  rocpts = [];
  prpts = [];
  nclus = size(simints,1);
  for topi = 1:(nclus^2-nclus) 
    rocpts(topi,1) = 1 - sum(~sortsimints(1:topi))/sum(~sortsimints); %% 1-FPR
    rocpts(topi,2) = sum(sortsimints(1:topi))/sum(sortsimints); %% TPR

    prpts(topi,1) = sum(sortsimints(1:topi))/sum(sortsimints); %% TPR
    prpts(topi,2) = sum(sortsimints(1:topi))/topi; %% PPV
  end

  rocpts = [ 1 0; rocpts; 0 1 ];
  prpts = [ 0 1; prpts; 1 0 ];

  %% figure; plot(rocpts(:,1),rocpts(:,2)); axis([0, 1, 0, 1],'manual');
  %% figure; plot(prpts(:,1),prpts(:,2)); axis([0, 1, 0, 1],'manual');

  auroc = -trapz(rocpts(:,1)',rocpts(:,2)');
  aupr = trapz(prpts(:,1)',prpts(:,2)');

end