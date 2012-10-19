function[ llh ] = geneInterLikelihood( D, par, priors, info )

expect = getExpectations(par);

%% intQlogP

mumat = par.mumean(info.time+1,:);
muSqExp = par.mumean.^2 + 1./(par.muprec);
muSqMat = muSqExp(info.time+1,:);
deltamat = repmat(par.deltamean',1,info.G);
deltaSqMat = repmat( par.deltamean'.^2 + 1./par.deltaprec', 1, info.G );

% sum(sum(D-mumat))
% sum(sum(D.^2-muSqMat))

Fprecmat = repmat(0, size(par.Fmean) );
for i = 1:size(par.Fmean,2)
    Fprecmat(:,i) = diag(par.Fprec{i});
end
expSqmshipF = (1./(Fprecmat')*par.mship) + (par.Fmean'*par.mship).^2;

partD1 = numel(D)/2 * ( -log(2*pi) + gamExpLogx(par.Lshape,par.Lscale) );
partD2 = sum(sum( -0.5*expect.L * (D.^2+muSqMat+deltaSqMat - 2*D.*deltamat - 2*D.*mumat + 2*mumat.*deltamat) ));

% blah1 = [ sum(sum(D.^2)) sum(sum(muSqMat)) sum(sum(deltaSqMat)) sum(sum(-2*D.*deltamat)) sum(sum(-2*D.*mumat)) sum(sum(2*mumat.*deltamat)) ];
% blahsum = -0.5*expect.L * sum(blah1)
% partD2 = blahsum;
% size(D.^2)
% size(muSqMat) 
% size(deltaSqMat)
% size(-2*D.*deltamat)
% size(-2*D.*mumat)
% size(2*mumat.*deltamat)

partMU1 = size(par.mumean,1) * sum( (-0.5*log(2*pi)+0.5*gamExpLogx(par.gamshape,par.gamscale)) );
partMU2 = sum(sum( -0.5*repmat(expect.gam,size(par.mumean(info.tuniq+1,:),1),1).*(muSqExp(info.tuniq+1,:) + expSqmshipF(info.tuniq+1,:) -2*par.mumean(info.tuniq+1,:).*(par.Fmean(:,info.tuniq+1)'*par.mship)) ));

%blah3 = sum(sum(muSqExp + expSqmshipF -2*par.mumean.*(par.Fmean'*par.mship)));
%-(sum(sum(muSqExp)) + sum(sum(expSqmshipF)) -2*sum(sum(par.mumean.*(par.Fmean'*par.mship))))

partDELTA1 = sum( -0.5*log(2*pi)+0.5*log(priors.deltaprec) );
partDELTA2 = sum( -0.5*priors.deltaprec.*(par.deltamean.^2+1./par.deltaprec - 2*par.deltamean.*priors.deltamean + priors.deltamean.^2) );

tempsum = repmat(0,1,info.nclus);
for i = 1:info.nclus
    temp = repmat(0,1,info.Tmax);
    for t = 1:info.Tmax
        tind = t+1;
        expFsq = par.Fmean(i,tind)^2 + 1/par.Fprec{tind}(i,i);
        
        % expSFsq = ( par.Smean(i,:)*par.Fmean(:,tind-1) )^2 + 1./(abs(par.Smean(i,:))*diag(par.Fprec{tind-1}));
        % expSFsq = (par.Smean(i,:)*par.Fmean(:,tind-1))^2 + ( ( par.Smean(i,:).^2+1./diag(par.Sprec{i})' ) * ( par.Fmean(:,tind-1).^2+1./(diag(par.Fprec{tind-1})) ) - par.Smean(i,:).^2*par.Fmean(:,tind-1).^2 );
        
        expSsq = par.Smean(i,:).^2 + 1./diag(par.Sprec{i})';
        expFminus1sq = par.Fmean(:,tind-1).^2 + 1./diag(par.Fprec{tind-1});
        SFelementwise = par.Smean(i,:).*par.Fmean(:,tind-1)';
        expSFsq = expSsq*expFminus1sq + sum(sum(SFelementwise'*SFelementwise)) - trace(SFelementwise'*SFelementwise);
         
        expSconstSq = par.Sconstmean(i)^2 + 1./par.Sconstprec(i);
        expFtSF = par.Fmean(i,tind) * par.Smean(i,:)*par.Fmean(:,tind-1);
        expFtSconst = par.Fmean(i,tind) * par.Sconstmean(i);
        expSFSconst = par.Smean(i,:)*par.Fmean(:,tind-1) * par.Sconstmean(i);
        temp(t) = expFsq +expSFsq +expSconstSq -2*expFtSF -2*expFtSconst +2*expSFSconst;
    end
    tempsum(i) = expect.sig(i,i) * sum(temp);
end

partF1prior1 = sum( -0.5*log(2*pi)+0.5*log(diag(priors.Fprec{1})) );
partF1prior2 = sum( -0.5*diag(priors.Fprec{1}).*( par.Fmean(:,1).^2+1./diag(par.Fprec{1})-2*par.Fmean(:,1).*priors.Fmean(:,1)+priors.Fmean(:,1).^2 ) );

partF2 = (info.Tmax-1)*sum( -0.5*log(2*pi) + 0.5*gamExpLogx(par.SIGshape,par.SIGscale) ) + sum( -0.5*tempsum );
partL = priors.Lshape*log(priors.Lscale) - gammaln(priors.Lshape) + (priors.Lshape-1)*gamExpLogx(par.Lshape,par.Lscale) - priors.Lscale*expect.L;
partGAM = sum( priors.gamshape.*log(priors.gamscale) - gammaln(priors.gamshape) + (priors.gamshape-1).*gamExpLogx(par.gamshape,par.gamscale) - priors.gamscale.*expect.gam );
partSIG = sum( priors.SIGshape.*log(priors.SIGscale) - gammaln(priors.SIGshape) + (priors.SIGshape-1).*gamExpLogx(par.SIGshape,par.SIGscale) - priors.SIGscale.*diag(expect.sig)' );

temp = repmat(0,1,info.nclus);
for j = 1:info.nclus
    tmp = -info.nclus/2*log(2*pi) + 0.5*logdet(priors.Sprec{j});
    temp(j) =  tmp - 0.5*( par.Smean(j,:)*par.Smean(j,:)'+sum(1./diag(par.Sprec{j})) -2*par.Smean(j,:)*priors.Smean(j,:)' + priors.Smean(j,:)*priors.Smean(j,:)' );
end
partS = sum(temp);

 QPparts = [ partD1 partD2 partMU1 partMU2 partDELTA1 partDELTA2 partF1prior1 partF1prior2 partF2 partL partGAM partSIG partS ...
%     sum(sum(D.^2)) ...
%     sum(sum(muSqMat)) ...
%     sum(sum(deltaSqMat)) ...
%     sum(sum(-2*D.*deltamat)) ...
%     sum(sum(-2*D.*mumat)) ...
%     sum(sum(2*mumat.*deltamat)) ...
%     expect.L ...
    ];
intQlogP = partD1+partD2+partMU1+partMU2+partDELTA1+partDELTA2+partF1prior1+partF1prior2+partF2+partL+partGAM+partSIG+partS;

%% int QlogQ (entropy of posteriors)

partMU = sum(sum( -gaussEnt(par.muprec) ));
partDELTA = sum( -gaussEnt(par.deltaprec) );

temp = repmat(0,1,info.Tmax+1);
for tind = 1:(info.Tmax+1)
    %par.Fprec{tind}
    temp(tind) = -multiGaussEnt(par.Fprec{tind},info.nclus);
end
partF = sum(temp);

partL = -gamEnt(par.Lshape,par.Lscale);

partGAM = sum( -gamEnt(par.gamshape,par.gamscale) );
partSIG = sum( -gamEnt(par.SIGshape,par.SIGscale) );

temp = repmat(0,1,info.nclus);
for j = 1:info.nclus
    temp(j) = -multiGaussEnt(par.Sprec{j},info.nclus);
end
partS = sum(temp);

QQparts = [ partMU partDELTA partF partL partGAM partSIG partS ];
intQlogQ = sum(QQparts); %partMU + partDELTA + partF + partL + partGAM + partSIG + partS;

%% put them together

%intQlogP
%intQlogQ
llh = intQlogP - intQlogQ;
%llh = [ intQlogP-intQlogQ intQlogP intQlogQ ];
%llh = partMU2;
%123
%partD2

if real(llh) ~= llh
    par.gamshape
    par.gamscale
    [ intQlogP-intQlogQ QPparts 0/0 QQparts ]
    error('llh not real')
end

end