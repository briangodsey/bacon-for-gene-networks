function[ par ] = updateMU( t, D, par, priors, info )

expect = getExpectations( par );

%% all at once

% par.muprec = repmat(expect.gam,size(par.mumean,1),1) + repmat(expect.L,size(par.mumean));
% 
% temp1 = repmat(expect.gam,size(par.mumean,1),1) .* (par.Fmean'*par.mship);
% temp2 = expect.L * ( D - repmat(par.deltamean',1,info.G) );
% 
% par.mumean = (temp1 + temp2) ./ par.muprec;


%% one t at a time
tnow = info.time == t;
tind = t+1;

%% testing
% mumat = par.mumean(info.time+1,:);
% muSqExp = par.mumean.^2; % + 1./(par.muprec);
% muSqMat = muSqExp(info.time+1,:);
% deltamat = repmat(par.deltamean',1,info.G);
% %deltaSqMat = repmat( par.deltamean'.^2 + 1./par.deltaprec', 1, info.G );
% deltaSqMat = repmat( par.deltamean'.^2, 1, info.G );
% closeToD1 = sum(sum(D.^2+muSqMat+deltaSqMat - 2*D.*deltamat - 2*D.*mumat + 2*mumat.*deltamat));
% absToD1 = sum(sum( (D-mumat-deltamat).^2 ));
% 
% %diffCloseAbs1 = closeToD1 - absToD1
% 
% partD2bef = sum(sum( -0.5*expect.L * (D.^2+muSqMat+deltaSqMat - 2*D.*deltamat - 2*D.*mumat + 2*mumat.*deltamat) ));
% 
% 
% Fprecmat = repmat(0, size(par.Fmean) );
% for i = 1:size(par.Fmean,2)
%     Fprecmat(:,i) = diag(par.Fprec{i});
% end
% expSqmshipF = (par.Fmean(:,info.tuniq+1)'*par.mship).^2;% + (1./(Fprecmat(:,info.tuniq+1)')*par.mship);
% closeToF1 = sum(sum( muSqExp + expSqmshipF -2*par.mumean.*(par.Fmean'*par.mship) ));
% absToF1 = sum(sum( (par.Fmean'*par.mship-par.mumean).^2 ));
% 
% partMU2beforeupdate = sum(sum( -0.5*repmat(expect.gam,size(par.mumean,1),1).*(muSqExp + expSqmshipF -2*par.mumean.*(par.Fmean'*par.mship)) ));

%% actual update
par.muprec(tind,:) = expect.gam + info.N(tind) * expect.L;

temp1 = expect.gam .* ( par.Fmean(:,tind)' * par.mship );
temp2 = expect.L * sum( D(tnow,:) - repmat( par.deltamean(tnow)',1,info.G), 1 );

par.mumean(tind,:) = (temp1 + temp2) ./ par.muprec(tind,:);

%% more testing
% mumat = par.mumean(info.time+1,:);
% muSqExp = par.mumean.^2;% + 1./(par.muprec);
% muSqMat = muSqExp(info.time+1,:);
% deltamat = repmat(par.deltamean',1,info.G);
% %deltaSqMat = repmat( par.deltamean'.^2 + 1./par.deltaprec', 1, info.G );
% deltaSqMat = repmat( par.deltamean'.^2, 1, info.G );
% closeToD2 = sum(sum( D.^2+muSqMat+deltaSqMat - 2*D.*deltamat - 2*D.*mumat + 2*mumat.*deltamat ));
% absToD2 = sum(sum( (D-mumat-deltamat).^2 ));
% 
% partD2aft = sum(sum( -0.5*expect.L * (D.^2+muSqMat+deltaSqMat - 2*D.*deltamat - 2*D.*mumat + 2*mumat.*deltamat) ));
% 
% Fprecmat = repmat(0, size(par.Fmean) );
% for i = 1:size(par.Fmean,2)
%     Fprecmat(:,i) = diag(par.Fprec{i});
% end
% expSqmshipF = (par.Fmean(:,info.tuniq+1)'*par.mship).^2;% + (1./(Fprecmat(:,info.tuniq+1)')*par.mship);
% closeToF2 = sum(sum( muSqExp + expSqmshipF -2*par.mumean.*(par.Fmean'*par.mship) ));
% absToF2 = sum(sum( (par.Fmean'*par.mship-par.mumean).^2 ));
% 
% partMU2afterupdate = sum(sum( -0.5*repmat(expect.gam,size(par.mumean,1),1).*(muSqExp + expSqmshipF -2*par.mumean.*(par.Fmean'*par.mship)) ));
% 
% % diffD2 = partD2aft-partD2bef
% % diffMU2 = partMU2afterupdate-partMU2beforeupdate
% % moveToD = closeToD2-closeToD1
% % moveToF = closeToF2-closeToF1
% % moveToDabs = absToD2-absToD1
% % moveToFabs = absToF2-absToF1

end