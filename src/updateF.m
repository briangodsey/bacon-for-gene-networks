function[ par ] = updateF( t, exper, D, par, priors, infor, Scovarray, temp )

  expect = getExpectations( par.exper{exper} );
  tind = t+1;
  tdatind = infor.Tdatind{exper}(infor.Tdat{exper}==t);

  %% update Fprec

  if t == 0
    %% if t is the first time point in the model
    ppart1 = priors.exper{exper}.Fprec1;
  else
    ppart1 = expect.SIG;
  end

  if t == infor.Tmax{exper}
    %% if t is the last time point in the model
    ppart2 = 0;
  else
    
    %%for k = 1:infor.nclus
    %%  Scovmat{k} = pinv(par.Sprec{k});
    %%end

    %% %% I put this outside the function because it was unnecessarily
    %% %% repeated for each t
    %% Scovarray = zeros(infor.nclus,infor.nclus,infor.nclus);
    %% for k = 1:infor.nclus
    %%   Scovarray(k,:,:) = pinv(par.Sprec{k});
    %% end

    %% %% middle, somewhat slow version
    %% temp2 = zeros(infor.nclus);
    %% for i = 1:infor.nclus
    %%   for j = 1:infor.nclus
    %% 	temp2(i,j) = diag(expect.SIG)'*Scovarray(:,i,j);
    %%   end
    %% end


    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PUT OUTSIDE THE FUNCTION
    %% because it was most of the compute time
    %% %% newest, fast version
    %% temp = zeros(infor.nclus);
    %% for i = 1:infor.nclus
    %%   temp(i,:) = diag(expect.SIG)'*reshape(Scovarray(:,i,:), ...
    %% 					    [infor.nclus infor.nclus 1]);
    %% end

    %% %% old, really slow version
    %%temp = zeros(infor.nclus);
    %%for i = 1:infor.nclus
    %%  for j = 1:infor.nclus
    %%	varmat = zeros(infor.nclus);
    %%	for k = 1:infor.nclus
    %%	  varmat(k,k) = Scovarray{k}(i,j);
    %%	end
    %%	temp(i,j) = trace(expect.SIG*varmat);
    %%  end
    %%end

    %%sum(sum(temp-temp2))

    ppart2 = par.Smean'*expect.SIG*par.Smean + temp;
  end

  if infor.N{exper}(tind)>0
    %% if the time point t has data
    ppart3 = diag(sum(par.mship,2) .* expect.gam(:,tdatind));
  else
    ppart3 = 0;
  end
				%4356
				%ppart1
				%ppart2
				%ppart3

  par.exper{exper}.Fprec{tind} = ppart1 + ppart2 + ppart3;

  if( sum(diag(par.exper{exper}.Fprec{tind})<0)>0 )
    %%     par.exper{exper}.Fprec{tind} = par.exper{exper}.Fprec{tind} ...
    %% 	+ 2*abs(min(diag(par.exper{exper}.Fprec{tind})));
    %% ppart1
    %% ppart2
    %% ppart3
    %% priors.exper{exper}
    %% par.exper{exper}.Fprec{tind}
    %% 				%par.Sprec
    error('diag(par.exper{exper}.Fprec{tind}) has a negative element')
  end


				%tind
				%logdet(par.Fprec{tind})
				% if(prod(diag(par.Fprec{tind}))==0)
				%     ppart1
				%     error('Determinant is zero.');
				% end

  %% update Fmean (uses current Fprec just calculated)
  if t == 0
    %% if t is the first time point in the model
    mpart1 = priors.exper{exper}.Fprec1 * priors.exper{exper}.Fmean1;
  else
    mpart1 = expect.SIG * (par.Smean * par.exper{exper}.Fmean(:,tind-1) ...
			   + par.Sconstmean);
  end

  if t == infor.Tmax{exper}
    %% if t is the last time point in the model
    mpart2 = 0;
  else
    mpart2 = par.Smean' * expect.SIG * (par.exper{exper}.Fmean(:,tind+1) ...
					- par.Sconstmean);
  end

  if infor.N{exper}(tind)>0
    %% if the time point t has data

    %% %% oldest, slow version
    %% mpart32 = zeros(infor.nclus,1);
    %% for g = 1:infor.G
    %%   for i=1:infor.nclus
    %% 	mpart32(i,1) = mpart32(i,1) + par.mship(i,g) ...
    %% 	    * ( expect.gam(i,tdatind) * par.exper{exper}.mumean(tdatind,g) );
    %%   end
    %% end
    %% 
    %% %% middle, not too fast version
    %% mpart33 = zeros(infor.nclus,1);
    %% for i=1:infor.nclus
    %%   mpart33(i,1) = par.mship(i,:) ...
    %% 	  * expect.gam(i,tdatind) * par.exper{exper}.mumean(tdatind,:)' ;
    %% end

    %% blazing fast version
    mpart3 = expect.gam(:,tdatind) .* ...
	(par.mship * par.exper{exper}.mumean(tdatind,:)');

    %% sum(sum(mpart3-mpart32))
    %% sum(sum(mpart33-mpart32))
    
  else
    mpart3 = 0;
  end

  %% t
  %% [ ppart1\mpart1-par.mship*par.mumean(tdatind,:)'/3 ppart1(:,1:5) ]
  %% [ ppart2\mpart2-par.mship*par.mumean(tdatind,:)'/3 ppart2(:,1:5) ]
  %% [ ppart3\mpart3-par.mship*par.mumean(tdatind,:)'/3 ppart3(:,1:5) ]
  %% par.Fprec{tind} 
  %% (mpart1+mpart2+mpart3)

  warning off
  par.exper{exper}.Fmean(:,tind) = par.exper{exper}.Fprec{tind} \ ...
      (mpart1+mpart2+mpart3);
  warning on

  %%[ par.exper{exper}.Fmean(:,tind) par.mship * par.exper{exper}.mumean(tdatind,:)' ]

end
