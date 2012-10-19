function[ par ] = updateSIG( exper, D, par, priors, infor )

  ## EFFICIENCY NOTE:
  %% The execution time for this function is dominated by 
  %% trace(expStS / Fprec) below and a little from invFprec up top, which 
  %% I think would be fairly difficult to make faster

  par = updateSIGscaleHYP( exper, D, par, priors, infor );
  expect = getExpectations( par.exper{exper} );

  par.exper{exper}.SIGshape = priors.exper{exper}.SIGshape ...
      + repmat(infor.Tmax{exper},1,infor.nclus)/2;


  %% prec-calculation for below
  invFprec = {};
  for tind = 2:(infor.Tmax{exper}+1)
    invFprec{tind} = pinv(par.exper{exper}.Fprec{tind});
  end

  %% %%%%%%%%%%%%%
  %% messy algorithm; can I make this faster?
  tempsum = repmat(0,1,infor.nclus);
  for i = 1:infor.nclus

    %% %% this is the old, slower version
    %% temp2 = repmat(0,6,infor.Tmax{exper}+1);
    %% for tind = 2:(infor.Tmax{exper}+1)
    %%   %%tind = t+1;
    %%   %%expFsq = par.exper{exper}.Fmean(i,tind)^2 + 1/par.exper{exper}.Fprec{tind}(i,i);
    %%   expFsq = par.exper{exper}.Fmean(i,tind)^2 + pinv(par.exper{exper}.Fprec{tind})(i,i);
    %% 
    %%   %% calculate the expectation of *SF(t-1)*
    %%   Fvec = par.exper{exper}.Fmean(:,tind-1);
    %%   Fprec = par.exper{exper}.Fprec{tind-1};
    %%   Svec = par.Smean(i,:);
    %%   Svec(Svec==0) = 0; % I don't know why this helps, but otherwise there could be Inf entries in Svec'*Svec
    %%   Sprec = par.Sprec{i};
    %%   expSF = Svec*Fvec;
    %%   %% simplecovs = (Fvec*Fvec' + pinv(Fprec)) ...
    %%   %%     .* (Svec'*Svec + pinv(Sprec)) ...
    %%   %%     - (Fvec.*Svec') * (Fvec'.*Svec);
    %%   %% varSF = sum(sum(simplecovs));
    %%   %% expSFsq = expSF^2 + varSF;
    %%   expStS = Svec'*Svec + pinv(Sprec);
    %%   expSFsq = Fvec'*expStS*Fvec + trace(expStS / Fprec);
    %% 
    %%   expSconstSq = par.Sconstmean(i)^2 + 1./par.Sconstprec(i);
    %%   expFtSF = par.exper{exper}.Fmean(i,tind) * par.Smean(i,:) * par.exper{exper}.Fmean(:,tind-1);
    %%   expFtSconst = par.exper{exper}.Fmean(i,tind) * par.Sconstmean(i);
    %%   expSFSconst = par.Smean(i,:)*par.exper{exper}.Fmean(:,tind-1) * par.Sconstmean(i);
    %%   
    %%   %% save the elements of the sum into a column for each t; adding up
    %%   %% similar numbers first helps avoid rounding errors; we later sum
    %%   %% the columns first and then the rows
    %%   temp2(:,tind) = [ expFsq; expSFsq; expSconstSq; -2*expFtSF; -2*expFtSconst; 2*expSFSconst ];
    %% end


    %% %%%%%%%%%%%%%
    %% newer, faster version
    Svec = par.Smean(i,:);
    Svec(Svec==0) = 0; % I don't know why this helps, but otherwise there could be Inf entries in Svec'*Svec
    Sprec = par.Sprec{i};
    expStS = Svec'*Svec + pinv(Sprec);
    expSconstSq = par.Sconstmean(i)^2 + 1./par.Sconstprec(i);

    temp = repmat(0,6,infor.Tmax{exper}+1);
    for tind = 2:(infor.Tmax{exper}+1)
      %%tind = t+1;
      %%expFsq = par.exper{exper}.Fmean(i,tind)^2 + 1/par.exper{exper}.Fprec{tind}(i,i);
      expFsq = par.exper{exper}.Fmean(i,tind)^2 + invFprec{tind}(i,i);
      %% calculate the expectation of *SF(t-1)*
      Fvec = par.exper{exper}.Fmean(:,tind-1);
      Fprec = par.exper{exper}.Fprec{tind-1};

      expSF = Svec*Fvec;
      expSFsq = Fvec'*expStS*Fvec + trace(expStS / Fprec);

      expFtSF = par.exper{exper}.Fmean(i,tind) * par.Smean(i,:) * par.exper{exper}.Fmean(:,tind-1);
      expFtSconst = par.exper{exper}.Fmean(i,tind) * par.Sconstmean(i);
      expSFSconst = par.Smean(i,:)*par.exper{exper}.Fmean(:,tind-1) * par.Sconstmean(i);
      
      %% save the elements of the sum into a column for each t; adding up
      %% similar numbers first helps avoid rounding errors; we later sum
      %% the columns first and then the rows
      temp(:,tind) = [ expFsq; expSFsq; expSconstSq; -2*expFtSF; -2*expFtSconst; 2*expSFSconst ];
    end

    
    tempsum(i) = sum(sum(temp,1));

    if(tempsum(i)<0)
      temp
      sum(temp,1)
      tempsum(i) = 1e-10
    end
  end

  par.exper{exper}.SIGscale = expect.SIGscaleHYP + 0.5 * tempsum;


  if( sum(sum(expect.SIG<0)) > 0 ) 
    %%expect = getExpectations( par.exper{exper} );
    %%par.Smean
    %%par.Sprec
    %%par.exper{exper}.Fmean
    %%par.exper{exper}.Fprec
    %% expect.SIG
    %% expect.SIGscaleHYP
    error('something in expect.SIG is negative')
  end


end