function[ par priors ] = updateVB( D, par, priors, infor, iter, doclustering )


  %%llhmat = zeros(size(D,2),5)
  %%llh.orig = geneInterLikelihood( D, par, priors, infor );
  
  %% llhmat2 = geneInterLikelihood( D, par, priors, infor );
  %% ab = llhmat2-llhmat1;
  %% llhmat1 = llhmat2;

  if mod(iter,5)==1
    if doclustering
      par = getMembership( par, priors, infor, 1, D );
    end
  end


  par = updateS( D, par, priors, infor );

  %% %%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% pre-calculated for updateF()
  %% it used to be inside the function but it was nearly 100%
  %% of the compute time, computed for each t. Duh.
  Scovarray = zeros(infor.nclus,infor.nclus,infor.nclus);
  for k = 1:infor.nclus
    Scovarray(k,:,:) = pinv(par.Sprec{k});
  end
  
  %% %%%%%%%%%%%%%%%%%%%%%
  %% loop over experiments
  for exper = 1:size(D,2)

    par = forwardPass( exper, D, par, priors, infor, Scovarray );
    par = updateGAM( exper, D, par, priors, infor );
    par = updateL( exper, D, par, priors, infor );

    
    if mod(iter,9)==1
      par = updateSIG( exper, D, par, priors, infor );
    end
    
    par = backwardPass( exper, D, par, priors, infor, Scovarray );
    par = updateGAM( exper, D, par, priors, infor );
    par = updateL( exper, D, par, priors, infor );
    %% par = updateSIG( exper, D, par, priors, infor );
    
  end

  %% ab
  %%llhmat-repmat(llhmat(1,:),size(llhmat,1),1)
  
  if iter > 150 && mod(iter,17)==1
    priors = updateHyperpriors( D, par, priors, infor );
  end

end
