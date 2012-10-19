function [ par ] = setpar( D, priors, info )

smallnum = 1e-4;


%% par.exper{exper}.SmeanMiddleGuymean = 2 * rand(info.nclus) - 1;
%% par.exper{exper}.SmeanMiddleGuyprec = maxPrec * rand(info.nclus);
%% par.exper{exper}.SmeanHYPmean = 2 * rand(1) - 1;
%% par.exper{exper}.SmeanHYPprec = maxPrec * rand(1);
%% par.exper{exper}.SmeanHYPdoublePRECshape = maxPrec * rand(1);
%% par.exper{exper}.SmeanHYPdoublePRECscale = maxPrec * rand(1);
%% par.exper{exper}.SprecHYPshape = maxPrec * rand(size(priors.SprecHYPshape));
%% par.exper{exper}.SprecHYPscale = maxPrec * rand(size(priors.SprecHYPscale));
par.Sconstmean = zeros(info.nclus,1);
par.Sconstprec = ones(info.nclus,1); %)maxPrec * rand(info.nclus,1);

%%par.Smean = eye(info.nclus); %% .* (rand(info.nclus) - 0.5);
par.Smean = ones(info.nclus) .* (rand(info.nclus) - 0.5);
for i = 1:info.nclus
  par.Sprec{i} = eye(info.nclus);
end



%%par.mship = rand(size(priors.mship));
[ s reorderedgenes ] = sort( rand(1,info.G) );
[ s reorderedclusts ] = sort( rand(1,info.nclus) );
remainingclusts = ceil( rand(1,info.G-info.nclus)*info.nclus );
%% par.mship = zeros(info.nclus,info.G);
ind = [ [ reorderedclusts' ; remainingclusts' ] reorderedgenes' ];
%%ind(:,3) = 1;
par.mship = full(sparse(ind(:,1),ind(:,2),1));


for exper = 1:size(D,2)
  dat = D{exper};

  mindat = min(min(dat));
  maxdat = max(max(dat));
  rangedat = maxdat-mindat;
  %%maxPrec = 10/rangedat;

  par.exper{exper}.mumean = rangedat * rand(length(info.Tdat{exper}),info.G) + mindat;
  par.exper{exper}.muprec = repmat(smallnum,length(info.Tdat{exper}),info.G); %maxPrec * rand(info.Tmax+1,info.G);

  par.exper{exper}.Lshape = priors.exper{exper}.Lshape; %maxPrec * rand(size(priors.Lshape));
  par.exper{exper}.Lscale = priors.exper{exper}.Lscale; %maxPrec * rand(size(priors.Lscale));


  %% the multivariate Gaussian cluster centers
  %% for Fmean, add T+1 for convenience; it stays zero
  par.exper{exper}.Fmean = rangedat * (rand(info.nclus,info.Tmax{exper}+1)) + mindat;
  par.exper{exper}.Fprec = {};
  for t = 1:info.Tmax{exper}+1
    par.exper{exper}.Fprec{t} = priors.exper{exper}.Fprec1;
  end

  %% the parameters of the precisions of the multivariate Gaussian clusters
  par.exper{exper}.gamshape = repmat( smallnum, info.nclus, length(info.Tdat{exper}) ); 
  par.exper{exper}.gamscale = ones(size(par.exper{exper}.gamshape));

  %% the hyperparameters of the beta/scale on the gamma dist of par.exper{exper}.gam
  par.exper{exper}.gamscaleHYPshape = smallnum;
  par.exper{exper}.gamscaleHYPscale = smallnum;

  par.exper{exper}.SIGscaleHYPshape = smallnum;
  par.exper{exper}.SIGscaleHYPscale = smallnum;

  par.exper{exper}.SIGshape = repmat(smallnum^2,1,info.nclus);
  par.exper{exper}.SIGscale = repmat(1,1,info.nclus);

end

end