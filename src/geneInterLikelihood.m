function[ llh ] = geneInterLikelihood( D, par, priors, infor )

  %% EFFICIENCY NOTE:
  %% this function is nearly 100% trace(expSS{i}/Fprec) by time
  %% which would be hard to make faster

  %% intQlogP

  partSconst = sum( -0.5*log(2*pi) + 0.5*log(priors.Sconstprec) ...
		   - 0.5*priors.Sconstprec ...
		   .* ( par.Sconstmean.^2 + 1./par.Sconstprec ...
		       - 2*par.Sconstmean.*priors.Sconstmean ...
		       + priors.Sconstmean.^2) );
  temp = repmat(0,1,infor.nclus);
  for j = 1:infor.nclus
    temp(j) = -infor.nclus/2*log(2*pi) +0.5*logdet(priors.Sprec{j}) ...
	- 0.5*( (par.Smean(j,:)-priors.Smean(j,:)) ...
	       * priors.Sprec{j} * (par.Smean(j,:)-priors.Smean(j,:))' ...
	       + trace( priors.Sprec{j} / par.Sprec{j} ) );
  end
  partS = partSconst + sum(temp);
  
  temp = par.mship .* log(priors.mship);
  temp(~priors.mship) = 0;
  partmship = sum(sum( temp ));
  %%partmship = sum(sum( par.mship .* log(priors.mship) ))


  %% %%%%%
  %% pre-calculation for below
  %%invSprec = {};
  expSS = {};
  for i = 1:infor.nclus
    invSprec = pinv(par.Sprec{i});
    expSS{i} = par.Smean(i,:)'*par.Smean(i,:) + invSprec;
  end
  

  for exper = 1:size(D,2)

    expect = getExpectations(par.exper{exper});

    mumat = par.exper{exper}.mumean(infor.timeD{exper},:);

    muSqExp = par.exper{exper}.mumean.^2 + 1./(par.exper{exper}.muprec);
    muSqMat = muSqExp(infor.timeD{exper},:);

    partD =  sum(sum(( -0.5*log(2*pi) + 0.5*gamExpLogx(par.exper{exper}.Lshape,par.exper{exper}.Lscale) ) ...
		     -0.5*expect.L * (D{exper}.^2 + muSqMat - 2*D{exper}.*mumat) ));

    Fmeans = par.exper{exper}.Fmean(:,infor.Tdat{exper}+1);
    Fprecmat = repmat(0, size(Fmeans) );
    for i = 1:size(Fprecmat,2)
      Fprecmat(:,i) = diag(par.exper{exper}.Fprec{infor.Tdat{exper}(i)+1});
    end
    expSqF = Fmeans.^2 + 1./Fprecmat;

    tempMU = zeros(1,infor.G);
    for g = 1:infor.G
      %% *temp* is a nclus x T matrix
      temp = par.mship(:,g)' ...
	  * ( -0.5*log(2*pi) ...
	     + 0.5*gamExpLogx(par.exper{exper}.gamshape, ...
			      par.exper{exper}.gamscale) ...
	     -  0.5*expect.gam ...
	     .* ( repmat(muSqExp(:,g)',infor.nclus,1) + expSqF ...
		 - 2*repmat(par.exper{exper}.mumean(:,g)',infor.nclus,1) ...
		 .*Fmeans ) );
      %% sum over all *t*
      tempMU(g) = sum(temp);
    end
    partMU = sum(tempMU);


    tempsumF = repmat(0,1,infor.Tmax{exper}+1);
    tempsummat = zeros(infor.Tmax{exper}+1,6);
    
    for tind = infor.Tall{exper}+1
      if tind==1
      	tempsumF(tind) = (par.exper{exper}.Fmean(:,tind)-priors.exper{exper}.Fmean1)' * priors.exper{exper}.Fprec1 ...
	    * (par.exper{exper}.Fmean(:,tind)-priors.exper{exper}.Fmean1) ...
	    + trace( priors.exper{exper}.Fprec1 / par.exper{exper}.Fprec{1} );
      else
	expFsq = par.exper{exper}.Fmean(:,tind)'*expect.SIG*par.exper{exper}.Fmean(:,tind) ...
	    + trace( expect.SIG / par.exper{exper}.Fprec{tind} );
	Fvec = par.exper{exper}.Fmean(:,tind-1);
	Fprec = par.exper{exper}.Fprec{tind-1};

	%% %% %%%%%%%%%%%%
	%% %% slow version
	%% tic
	%% for i = 1:infor.nclus
	%%   invSprec = pinv(par.Sprec{i});
	%%   for j = 1:infor.nclus
	%%     expSS = par.Smean(i,:)'*par.Smean(j,:);
	%%     if i == j
	%%       expSS = expSS + invSprec;
	%%     end
	%%     covSF(i,j) = Fvec'*expSS*Fvec + trace(expSS/Fprec) ...
	%% 	- par.Smean(i,:)*Fvec*par.Smean(j,:)*Fvec;
	%%   end
	%% end
	%% toc
	%% expSFsq1 = (par.Smean*Fvec)'*expect.SIG*(par.Smean*Fvec) ...
	%%     + trace( expect.SIG * covSF );

	%% %%%%%%%%%%%%%%
	%% faster version where we calculate only the diagonal
	%% because that's all that we need later
	covSF = zeros(infor.nclus);
	for i = 1:infor.nclus
	  %% invSprec = pinv(par.Sprec{i});
	  %% expSS = par.Smean(i,:)'*par.Smean(i,:) + invSprec;
	  covSF(i,i) = Fvec'*expSS{i}*Fvec + trace(expSS{i}/Fprec) ...
	      - par.Smean(i,:)*Fvec*par.Smean(i,:)*Fvec;
	end

	expSFsq = (par.Smean*Fvec)'*expect.SIG*(par.Smean*Fvec) ...
	    + trace( expect.SIG * covSF );

	%%abs(expSFsq-expSFsq1)

	expSconstSq = par.Sconstmean'*expect.SIG*par.Sconstmean ...
	    + trace(expect.SIG/diag(par.Sconstprec));

	expFtSF = par.exper{exper}.Fmean(:,tind)' * expect.SIG * par.Smean*par.exper{exper}.Fmean(:,tind-1);
	expFtSconst = par.exper{exper}.Fmean(:,tind)' * expect.SIG * par.Sconstmean;
	expSFSconst = (par.Smean*par.exper{exper}.Fmean(:,tind-1))' * expect.SIG * par.Sconstmean;

	tempsumF(tind) = expFsq +expSFsq +expSconstSq ...
	    -2*expFtSF -2*expFtSconst +2*expSFSconst;

	tempsummat(tind,:) = [ expFsq expSFsq expSconstSq expFtSF expFtSconst expSFSconst];
      end
    end

    tempF = zeros(1,infor.Tmax{exper}+1);
    tempF(1) = -infor.nclus/2*log(2*pi) ...
	+ 0.5*logdet(priors.exper{exper}.Fprec1) ...
	- 0.5*tempsumF(1);
    for tind = 2:(infor.Tmax{exper}+1)
      tempF(tind) = -infor.nclus/2*log(2*pi) ...
	  + 0.5*sum(gamExpLogx(par.exper{exper}.SIGshape,par.exper{exper}.SIGscale)) ...
	  - 0.5*tempsumF(tind);
    end
    partF = sum(tempF);

    partL = priors.exper{exper}.Lshape*log(priors.exper{exper}.Lscale) - gammaln(priors.exper{exper}.Lshape) ...
	+ (priors.exper{exper}.Lshape-1)*gamExpLogx(par.exper{exper}.Lshape,par.exper{exper}.Lscale) ...
	- priors.exper{exper}.Lscale*expect.L;
    
    partgamscaleHYP = priors.exper{exper}.gamscaleHYPshape*log(priors.exper{exper}.gamscaleHYPscale) ...
	- gammaln(priors.exper{exper}.gamscaleHYPshape)  ...
	+ (priors.exper{exper}.gamscaleHYPshape-1)*gamExpLogx(par.exper{exper}.gamscaleHYPshape,par.exper{exper}.gamscaleHYPscale) ...
	- priors.exper{exper}.gamscaleHYPscale*expect.gamscaleHYP;
    partGAM = partgamscaleHYP ...
	+ sum(sum( priors.exper{exper}.gamshape.*gamExpLogx(par.exper{exper}.gamscaleHYPshape,par.exper{exper}.gamscaleHYPscale) ...
		  - gammaln(priors.exper{exper}.gamshape) ...
		  + (priors.exper{exper}.gamshape-1).* gamExpLogx(par.exper{exper}.gamshape,par.exper{exper}.gamscale) ...
		  - expect.gamscaleHYP.*expect.gam ));

    partSIGscaleHYP = priors.exper{exper}.SIGscaleHYPshape*log(priors.exper{exper}.SIGscaleHYPscale) ...
	- gammaln(priors.exper{exper}.SIGscaleHYPshape)  ...
	+ (priors.exper{exper}.SIGscaleHYPshape-1)*gamExpLogx(par.exper{exper}.SIGscaleHYPshape,par.exper{exper}.SIGscaleHYPscale) ...
	- priors.exper{exper}.SIGscaleHYPscale*expect.SIGscaleHYP;
    partSIG = partSIGscaleHYP + ...
	sum( priors.exper{exper}.SIGshape.*gamExpLogx(par.exper{exper}.SIGscaleHYPshape,par.exper{exper}.SIGscaleHYPscale) ...
	    - gammaln(priors.exper{exper}.SIGshape) ...
	    + (priors.exper{exper}.SIGshape-1).*gamExpLogx(par.exper{exper}.SIGshape,par.exper{exper}.SIGscale) ...
	    - expect.SIGscaleHYP.*diag(expect.SIG)' );

    QPparts(exper) = sum( [ partD partMU partF partL partGAM partSIG ] );
    
  end
  
  intQlogP = sum(QPparts) + partS + partmship;
  
  

  %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% int QlogQ (entropy of posteriors)

  %% partSmeanMiddleGuy  = sum(sum( -gaussEnt(par.SmeanHYPprec) ));
  %% partSmeanHYPdoublePREC = -gamEnt(par.SmeanHYPdoublePRECshape, ...
  %%     par.SmeanHYPdoublePRECscale);
  %% partSmeanHYP = -gaussEnt(par.SmeanHYPprec);
  %% partSprecHYP = -gamEnt(par.SprecHYPshape,par.SprecHYPscale);
  %% partSconst = sum( -gaussEnt(par.Sconstprec) );
  %% temp = repmat(0,1,infor.nclus);
  %% for j = 1:infor.nclus
  %%     temp(j) = -multiGaussEnt(par.Sprec{j});
  %% end
  %% partS = partSmeanMiddleGuy + partSmeanHYPdoublePREC + partSmeanHYP ...
  %%     + partSprecHYP + partSconst + sum(temp);

  partSconst = sum( -gaussEnt(par.Sconstprec) );
  temp = repmat(0,1,infor.nclus);
  for j = 1:infor.nclus
    temp(j) = -multiGaussEnt(par.Sprec{j});
  end
  partS = partSconst + sum(temp);

  temp = par.mship .* log(par.mship);
  temp(~par.mship) = 0;
  partmship = sum(sum( temp ));

  for exper = 1:size(D,2)
    expect = getExpectations(par.exper{exper});

    partMU = sum(sum( -gaussEnt(par.exper{exper}.muprec) ));

    temp = repmat(0,1,infor.Tmax{exper}+1);
    for tind = infor.Tall{exper}+1
      %%par.exper{exper}.Fprec{tind}
      temp(tind) = -multiGaussEnt(par.exper{exper}.Fprec{tind});
    end
    partF = sum(temp);

    partL = -gamEnt(par.exper{exper}.Lshape,par.exper{exper}.Lscale);

    partGAMscaleHYP = -gamEnt(par.exper{exper}.gamscaleHYPshape,par.exper{exper}.gamscaleHYPscale);
    partGAM = partGAMscaleHYP + sum(sum( -gamEnt(par.exper{exper}.gamshape,par.exper{exper}.gamscale) ));

    partSIGscaleHYP = -gamEnt(par.exper{exper}.SIGscaleHYPshape,par.exper{exper}.SIGscaleHYPscale);
    partSIG = partSIGscaleHYP + sum( -gamEnt(par.exper{exper}.SIGshape,par.exper{exper}.SIGscale) );



    QQparts(exper) = sum( [ partMU partF partL partGAM partSIG ] );

  end

  intQlogQ = sum(QQparts) + partS + partmship;

  %% put them together

  %%QPparts
  %%QQparts
  %%[ intQlogP intQlogQ ]
  
  llh = intQlogP - intQlogQ;
  %%llh = [ intQlogP-intQlogQ QPparts(1) QPparts(2:7)-QQparts ];
  %%llh = [ intQlogP-intQlogQ QPparts(1) QPparts(2:7)-QQparts sum(tempsummat,1) ];
  
  %%sum(tempsummat,1)
  
				%llh = partMU2;
				%123
  
  %% if real(llh) ~= llh
  %%   par.exper{exper}.gamshape
  %%   par.exper{exper}.gamscale
  %%   par.exper{exper}.gamscaleHYPshape
  %%   par.exper{exper}.gamscaleHYPscale
  %%   par.exper{exper}.SIGshape
  %%   par.exper{exper}.SIGscale
  %%   [ intQlogP-intQlogQ QPparts 0/0 QQparts ]
  %%   error('llh not real')
  %% end

end