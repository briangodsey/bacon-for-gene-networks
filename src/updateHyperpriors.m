function[ priors ] = updateHyperpriors( D, par, priors, infor )


  
  %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  priors.Sconstprec = ones(infor.nclus,1) * ...
      infor.nclus / ( (par.Sconstmean-priors.Sconstmean)' ...
		    * (par.Sconstmean-priors.Sconstmean) ...
		    + trace(pinv(diag(par.Sconstprec))) );

  %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  tempsum = 0;
  for j = 1:infor.nclus
    tempsum = tempsum + ...
	infor.nclus / ( (par.Smean(j,:)-priors.Smean(j,:)) ...
		      * (par.Smean(j,:)-priors.Smean(j,:))' ...
		      + trace(pinv(par.Sprec{j})) );
  end

  tmaxsum = 0;
  for exper = 1:size(D,2)
    tmaxsum = tmaxsum + infor.Tmax{exper};
  end

  for j = 1:infor.nclus
    priors.Sprec{j} = eye(infor.nclus) * ...
	tempsum / (tmaxsum+1);
  end


  %% loop over experiments
  for exper = 1:size(D,2)

    expect = getExpectations(par.exper{exper});

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    priors.exper{exper}.Fprec1 = eye(infor.nclus) * ...
    	infor.nclus / ( (par.exper{exper}.Fmean(:,1)-priors.exper{exper}.Fmean1)' ...
    		      * (par.exper{exper}.Fmean(:,1)-priors.exper{exper}.Fmean1) ...
    		      + trace(pinv(par.exper{exper}.Fprec{1})) );

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    priors.exper{exper}.Lshape = fsolve(@(x) priorHYPshapederiv(x, priors.exper{exper}.Lscale, ...
						   par.exper{exper}.Lshape, ...
						   par.exper{exper}.Lscale),
			   priors.exper{exper}.Lshape);
    priors.exper{exper}.Lscale = priors.exper{exper}.Lshape/expect.L;

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% UPDATE *GAM* STUFF

    %%psiroot = fsolve(@psi,1.4616);
    %%psiroot = 1.461632143;

    priors.exper{exper}.gamscaleHYPshape = fsolve(@(x) priorHYPshapederiv(x, priors.exper{exper}.gamscaleHYPscale, ...
							     par.exper{exper}.gamscaleHYPshape, ...
							     par.exper{exper}.gamscaleHYPscale),
				     priors.exper{exper}.gamscaleHYPshape);
    priors.exper{exper}.gamscaleHYPscale = priors.exper{exper}.gamscaleHYPshape / expect.gamscaleHYP;

    priors.exper{exper}.gamshape = fsolve(@(x) priorshapederiv(x, par.exper{exper}.gamscaleHYPshape, ...
						  par.exper{exper}.gamscaleHYPscale, ...
						  par.exper{exper}.gamshape,par.exper{exper}.gamscale),
			     priors.exper{exper}.gamshape);

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% UPDATE *SIG* STUFF

    priors.exper{exper}.SIGscaleHYPshape = fsolve(@(x) 
						  priorHYPshapederiv(x, priors.exper{exper}.SIGscaleHYPscale, ...
								     par.exper{exper}.SIGscaleHYPshape, ...
								     par.exper{exper}.SIGscaleHYPscale),
				     priors.exper{exper}.SIGscaleHYPshape);
    priors.exper{exper}.SIGscaleHYPscale = priors.exper{exper}.SIGscaleHYPshape / expect.SIGscaleHYP;

    priors.exper{exper}.SIGshape = fsolve(@(x) priorshapederiv(x, par.exper{exper}.SIGscaleHYPshape, ...
						  par.exper{exper}.SIGscaleHYPscale, ...
						  par.exper{exper}.SIGshape,par.exper{exper}.SIGscale),
			     priors.exper{exper}.SIGshape);

    %% priors.exper{exper}.gamshape
    %% priors.exper{exper}.gamscaleHYPshape
    %% priors.exper{exper}.SIGshape
    %% priors.exper{exper}.SIGscaleHYPshape



  end

end