function[ par ] = getMembership( par, priors, infor, shuffle, D )

  warning("off","all")

  %% %% check to see if two clusters are the same
  %% for i = 2:infor.nclus
  %%   for j = 1:i-1
  %%     if( all( abs(par.Fmean(i,:)./par.Fmean(j,:)-1)<0.05 ) )
  %%        par.Fmean(i,:) = 2*(rand(1,size(par.Fmean,2))-0.5);
  %%     end
  %%   end
  %% end


  temp = {};
  %%temp2 = {};
  %% %%%%%%%%%%%%%%%%%%%%%
  %% loop over experiments
  for exper = 1:size(D,2)
    
    expect = getExpectations( par.exper{exper} );
    
    %% %% %%%%%%%%%%%
    %% %% old, slow version
    %% tic
    %% temp2{exper} = zeros(infor.nclus,infor.G);
    %% for g = 1:infor.G
    %%   mus = par.exper{exper}.mumean(:,g);
    %%   muprecs = diag(par.exper{exper}.muprec(:,g));
    %%   for i = 1:infor.nclus
    %% 	expgam = diag(expect.gam(i,:));
    %% 	Fmeans = par.exper{exper}.Fmean(i,infor.Tdat{exper}+1)';
    %% 	Fprecs = zeros(length(infor.Tdat{exper}));
    %% 	for tind = infor.Tdatind{exper}
    %% 	  Fprecs(tind,tind) = par.exper{exper}.Fprec{infor.Tdat{exper}(tind)+1}(i,i);
    %% 	end
    %% 
    %% 	temp2{exper}(i,g) = (-0.5*length(infor.Tdat{exper}))*log(2*pi) ...
    %% 	    + 0.5* sum(gamExpLogx(par.exper{exper}.gamshape(i,:), ...
    %% 				  par.exper{exper}.gamscale(i,:))) ...
    %% 	    - 0.5 * ( mus'*expgam*mus + trace(expgam/muprecs) ...
    %% 		     + Fmeans'*expgam*Fmeans + trace(expgam/Fprecs) ...
    %% 		     - 2*mus'*expgam*Fmeans );
    %%   end
    %% end

    %% %%%%%%%%%%%%%%
    %% newer, faster version
    temp{exper} = zeros(infor.nclus,infor.G);
    for i = 1:infor.nclus
      Fmeans = par.exper{exper}.Fmean(i,infor.Tdat{exper}+1)';
      
      Fprecs = zeros(length(infor.Tdat{exper}));
      for tind = infor.Tdatind{exper}
	Fprecs(tind,tind) = par.exper{exper}.Fprec{infor.Tdat{exper}(tind)+1}(i,i);
      end
     
      temp{exper}(i,:) = (-0.5*length(infor.Tdat{exper}))*log(2*pi) ...
	  + 0.5* sum(gamExpLogx(par.exper{exper}.gamshape(i,:), ...
				par.exper{exper}.gamscale(i,:))) ...
	  - 0.5 * ( diag( par.exper{exper}.mumean' ...
			 * diag(expect.gam(i,:)) * par.exper{exper}.mumean )' ...
		   + sum( diag(expect.gam(i,:)) * par.exper{exper}.muprec.^-1 ) ...
		   + Fmeans'* diag(expect.gam(i,:)) *Fmeans ...
		   + trace( diag(expect.gam(i,:)) / Fprecs ) ...
		   - 2 * Fmeans' * diag(expect.gam(i,:)) ...
		     * par.exper{exper}.mumean );
      
    end

    %%sum(sum(temp{exper}-temp2{exper}))

  end


  tempsum = zeros(size(temp{1}));
  for exper=1:size(temp,2)
    tempsum = tempsum + temp{exper};
  end

  logllh = log(priors.mship) + real(tempsum);

  %% %% %%%%%%%%%%%%%%%%
  %% %% older, slower version
  %% tic
  %% mship2 = zeros( infor.nclus, infor.G );
  %% %%exp(logllh)
  %% for g = 1:infor.G
  %%   llhrow = logllh(:,g);
  %%   grow = exp(llhrow-max(llhrow));
  %%   %% if all(grow==0) || sum(grow>0.99*max(grow))>1
  %%   %%   grow(logllh(:,g)==max(logllh(:,g))) = 1;
  %%   %% end
  %%   %% if sum( grow>0.99*max(grow) ) > 1
  %%   %%   ind = find( grow>0.95*max(grow) );
  %%   %%   grow(ind(mod(g,length(ind))+1)) = 2*max(grow);
  %%   %% end
  %%   mship2(:,g) = grow/sum(grow);
  %% end
  %% toc

  %%tic
  mshipUNNORM = exp( logllh - repmat(max(logllh),infor.nclus,1) );
  mship = mshipUNNORM ./ repmat(sum(mshipUNNORM),infor.nclus,1);
  %%toc

  %%sum(sum(abs(mship-mship2)))


  if shuffle
    %% if a cluster is empty
    clustmemb = sum(mship>0.1,2);
    emptyclusters = find( clustmemb==0  )';
    if( size(emptyclusters,2) ~= 0 )
      numemptyclusters = length(emptyclusters)
      for i = emptyclusters
	clustmemb = sum(mship,2);
	biggestclust = find(clustmemb>mean(clustmemb));
	bigclustermemb = find(max(mship(biggestclust,:),[],1)>0.1);
	maxlogllh = max(logllh);
	%%loneliestgene = find( maxlogllh == min(maxlogllh(clustmemb==max(clustmemb))) );
	loneliestgene = find(maxlogllh==min(maxlogllh(bigclustermemb)));
	mship(:,loneliestgene) = 0;
	mship(i,loneliestgene) = 1;
	maxlogllh(loneliestgene) = max(maxlogllh)+1;
	%%clustmemb = sum(mship,2)

	%% make the lonely gene comfortable in its new home
	%% loop over experiments

  %% %%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% pre-calculated for updateF()
  %% it used to be inside the function but it was nearly 100%
  %% of the compute time, computed for each t. Duh.
  Scovarray = zeros(infor.nclus,infor.nclus,infor.nclus);
  for k = 1:infor.nclus
    Scovarray(k,:,:) = pinv(par.Sprec{k});
  end

	for exper = 1:size(D,2)
	  par.exper{exper}.Fmean(i,:) = par.exper{exper}.mumean(:,loneliestgene)';
	  par.exper{exper}.gamshape(i,:) = max(max(par.exper{exper}.gamshape));
	  par.exper{exper}.gamscale(i,:) = min(min(par.exper{exper}.gamscale));
	  par = forwardPass( exper, D, par, priors, infor, Scovarray );
	  %%par = backwardPass( exper, D, par, priors, infor, Scovarray  );
	end
      end
    end
  end

  par.mship = mship;

end
