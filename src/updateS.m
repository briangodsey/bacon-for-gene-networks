function[ par ] = updateS( D, par, priors, infor )
  
  par = updateSconst( D, par, priors, infor );

  ## pre-calculations for temp1 below
  temp0 = {};
  for exper = 1:size(D,2)
    %% Finvmat = zeros(infor.nclus, infor.nclus, ...
    %% 		    infor.Tmax{exper}+1);
    %% for tsum = 2:(infor.Tmax{exper}+1)
    %%   Finvmat(:,:,tsum-1) = pinv(par.exper{exper}.Fprec{tsum-1});
    %% end
    temp0{exper} = zeros(infor.nclus);
    for tsum = 2:(infor.Tmax{exper}+1)
      temp0{exper} = temp0{exper} + ...
	  ( par.exper{exper}.Fmean(:,tsum-1) ...
	   * par.exper{exper}.Fmean(:,tsum-1)' ... 
	   + pinv(par.exper{exper}.Fprec{tsum-1}) );
	   %%+ Finvmat(:,:,tsum-1) );
    end
  end

  for i = 1:infor.nclus

    temp1 = {};
    sum2 = {};
    %% %%%%%%%%%%%%%%%%%%%%%
    %% loop over experiments
    for exper = 1:size(D,2)
      
      expect = getExpectations( par.exper{exper} );
      
      %% %%%%%%%%%%%%%%%%%%%%
      %% calculate the row of the precision for this cluster

      %% %% can I speed this up?
      %% temp12{exper} = zeros(infor.nclus);
      %% for tsum = 2:(infor.Tmax{exper}+1)
      %% 	temp12{exper} = temp12{exper} + expect.SIG(i,i) ...
      %% 	    * ( par.exper{exper}.Fmean(:,tsum-1) ...
      %% 	       * par.exper{exper}.Fmean(:,tsum-1)' ... 
      %% 	       + pinv(par.exper{exper}.Fprec{tsum-1}) );
      %% end
      %%
      %% yes, you idiot; see the pre-calculation above used here:
      temp1{exper} = expect.SIG(i,i) * temp0{exper};


      if( sum(diag(par.Sprec{i})<0)>0 )
	i
	expect.SIG
	sum(expect.SIG<0)
	par.Sprec{i}
				%%par.Fprec
	error('diag(par.Sprec{i}) has a negative element');
      end
      
      %% %%%%%%%%%%%%%%%%%%%%%%%%
      %% calculate the mean for this cluster
      
      %% %% old, slower version
      %% temp2 = repmat(0,infor.Tmax{exper}+1,infor.nclus);
      %% for tsum = 2:(infor.Tmax{exper}+1)
      %% 	temp2(tsum,:) = (par.exper{exper}.Fmean(i,tsum)-par.Sconstmean(i)) ...
      %% 	    * expect.SIG(i,i) * par.exper{exper}.Fmean(:,tsum-1)';
      %% end
      %% sum2{exper} = sum(temp2,1);


      %% newer, faster version
      sum2{exper} = expect.SIG(i,i) ...
	  * (par.exper{exper}.Fmean(i,2:end)-par.Sconstmean(i)) ...
	  *  par.exper{exper}.Fmean(:,(2:end)-1)';

      %%sum(sum(sum2{exper}-sum22{exper}))
      
    end

    %% add up values from different experiments
    temp1sum = zeros(size(temp1{1}));
    sum2sum = zeros(size(sum2{1}));
    for exper=1:size(temp1,2)
      temp1sum = temp1sum + temp1{exper};
      sum2sum = sum2sum + sum2{exper};
    end

    par.Sprec{i} = priors.Sprec{i} + temp1sum;
    par.Smean(i,:) =  ( priors.Smean(i,:)*priors.Sprec{i} + sum2sum ) ...
	/ par.Sprec{i} ;

  end

end