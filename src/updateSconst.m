function[ par ] = updateSconst( D, par, priors, info )

  for i = 1:info.nclus
    
    temp1 = {};
    sum2 = {};
    %% %%%%%%%%%%%%%%%%%%%%%
    %% loop over experiments
    for exper = 1:size(D,2)
      
      expect = getExpectations( par.exper{exper} );
      
      %% %%%%%%%%%%%%%%%%%
      %% calculate the row of the precision for this cluster
      
      %% %% old, slower version
      %% temp12{exper} = 0;
      %% for tsum = 2:(info.Tmax{exper}+1)
      %%   temp12{exper} = temp12{exper} + expect.SIG(i,i);
      %% end

      %% newer, faster version
      temp1{exper} =  info.Tmax{exper} * expect.SIG(i,i);

      %%sum(sum(temp1{exper}-temp12{exper}))


      %% %%%%%%%%%%%%%%%%%
      %% calculate the mean for this cluster

      %% %% older, slow version
      %% temp2 = repmat(0,info.Tmax{exper}+1,1);
      %% for tsum = 2:(info.Tmax{exper}+1)
      %%   temp2(tsum) = (par.exper{exper}.Fmean(i,tsum)-par.Smean(i,:) ...
      %% 		       *par.exper{exper}.Fmean(:,tsum-1)) * expect.SIG(i,i);
      %% end
      %% sum22{exper} = sum(temp2);

      %% newer, faster version
      sum2{exper} = expect.SIG(i,i) ...
	  * sum( par.exper{exper}.Fmean(i,2:end) ...
		- par.Smean(i,:) * par.exper{exper}.Fmean(:,(2:end)-1) );

      %% sum(sum(sum2{exper}-sum22{exper}))


      
    end

    %% add up values from different experiments
    temp1sum = zeros(size(temp1{1}));
    sum2sum = zeros(size(sum2{1}));
    for exper=1:size(temp1,2)
      temp1sum = temp1sum + temp1{exper};
      sum2sum = sum2sum + sum2{exper};
    end

    par.Sconstprec(i) = priors.Sconstprec(i) + temp1sum;
    par.Sconstmean(i) =  ( priors.Sconstmean(i)*priors.Sconstprec(i) + sum2sum ) ...
	/ par.Sconstprec(i);

  end

end