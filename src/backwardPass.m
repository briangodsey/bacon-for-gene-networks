function[ par ] = backwardPass( exper, D, par, priors, infor, Scovarray )

  expect = getExpectations( par.exper{exper} );

  %% a pre-compute for updateF()
  %% obvs, it doesn't depend on t
  temp = zeros(infor.nclus);
  for i = 1:infor.nclus
    temp(i,:) = diag(expect.SIG)'*reshape(Scovarray(:,i,:), ...
  					    [infor.nclus infor.nclus 1]);
  end

  for t = fliplr( infor.Tall{exper} )
    if( sum(infor.time{exper} == t) >= 1 )
      par =  updateMU( t, exper, D, par, priors, infor );
    end
    par =  updateF( t, exper, D, par, priors, infor, Scovarray, temp );
  end

end