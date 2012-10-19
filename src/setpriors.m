function [priors] = setpriors( D, info )

  smallnum = 1e-2;
  largenum = 1e2;

  %% priors.SmeanHYPmean = 0;
  %% priors.SmeanHYPprec = smallnum;
  %% priors.SmeanHYPdoublePRECshape = 1;
  %% priors.SmeanHYPdoublePRECscale = 1;
  %% priors.SprecHYPshape = 1;
  %% priors.SprecHYPscale = 1;

  priors.Smean = zeros(info.nclus);
  priors.Sprec = {};
  for i = 1:info.nclus
    priors.Sprec{i} = eye(info.nclus);
  end

  priors.Sconstmean = zeros(info.nclus,1);
  priors.Sconstprec = ones(info.nclus,1);

  priors.mship = repmat(1/info.nclus, info.nclus, info.G);


  %% loop over experiments
  for exper = 1:size(D,2)

    %%priors.exper{exper}.mumean = repmat(0,info.Tmax+1,info.G);
    %%priors.exper{exper}.muprec = repmat(smallnum,info.Tmax+1,info.G);

    priors.exper{exper}.Lshape = smallnum;
    priors.exper{exper}.Lscale = smallnum;


    %% the multivariate Gaussian cluster centers
    %% for Fmean, add T+1 for convenience; it stays zero
    priors.exper{exper}.Fmean1 = zeros(info.nclus, 1);
    priors.exper{exper}.Fprec1 = smallnum*eye(info.nclus);

    %% the parameters of the precisions of the multivariate Gaussian clusters
    %%priors.exper{exper}.gamshape = repmat( 1, info.nclus, length(info.Tdat) );
    priors.exper{exper}.gamshape = 1;

    %% the hyperparameters of the beta/scale on the gamma dist of par.gam
    priors.exper{exper}.gamscaleHYPshape = smallnum; %info.Tmax;
    priors.exper{exper}.gamscaleHYPscale = smallnum; %info.Tmax;

    %%priors.exper{exper}.SIGshape = repmat(smallnum^2,1,info.nclus);
    priors.exper{exper}.SIGshape = 1;

    priors.exper{exper}.SIGscaleHYPshape = smallnum;
    priors.exper{exper}.SIGscaleHYPscale = smallnum;

  end

end