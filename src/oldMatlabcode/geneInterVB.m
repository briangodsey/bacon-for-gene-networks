addpath('~psykacek/public/mlbsrc/slide_normalization')
addpath('~bgodsey/AmpliMethods')
addpath('~bgodsey/src/matlab')

% [H, rawD] = scan_fspma('/bi/home/bgodsey/camda2008/printDataMatrix.tsv');
% rawD( rawD<=0 ) = min(min( rawD(rawD>0) ))/2;
% D = log( rawD(:,1:1000) );
% % clusternumbers = 10;


% the model always starts with t=0. If you want a prior state on the first time
% point, start the data at t=1. If the data starts at t=0, there are simply
% fixed weakly informative priors on time point t=0. Starting with data at
% t=0 seems to work better
%timeORIG = [ 0:7 0:7 0:7 ];
% %timeORIG = [ 0 0.5 1.5 3 6 9 12 24 0 0.5 1.5 3 6 9 12 24 0 0.5 1.5 3 6 9 12 24 ];

exper = repmat( 1, 1, 24 );
names = { 'apoptosis' };

simdata;
timeORIG = repmat( 0:(size(preD,2)-1), 1, numreps );

figure; plot(preD'); title('Simulated Data');
%print -depsc 

info = setinfo(D,timeORIG,exper,names);

%% Does the model fitting for various cluster numbers
clusternumbers = 5;
hyperparvals = 0.00001; %(1:20:200) / ( max(max(D))-min(min(D)) );
llh = repmat( -10^20, size(hyperparvals,2), size(clusternumbers,2) );
for i = 1:size(clusternumbers,2)
    nclus = clusternumbers(i);
    nclus
    save nclus.out nclus -ASCII -tabs

    info.nclus = nclus;
    priors = setpriors( D, info );
    %par = priors;
    %expect = getExpectations( par );

    for j = 1:size(hyperparvals,2)
        hp = hyperparvals(j);
        for t = 1:info.Tmax+1
            priors.Fprec{t} = eye( info.nclus )*hp;
        end
        [ par llh(j,i) iter ] = vbmodelfit( D, priors, info );
    end
end

% figure; plot(hyperparvals',llh);
% legend(num2str(clusternumbers'));
% title('Log likelihoods for various hyperparameter values');
% 
% figure; plot(hyperparvals',llh(:,1:3));
% legend(num2str(clusternumbers(1:3)'));
% title('Subset of log likelihoods for various hyperparameter values');

sortllh = sort(llh,'descend');
bestllh = [];
if(size(clusternumbers,2)>1)
    for i = 1:min(10,size(clusternumbers,2))
        matches = find(llh == sortllh(i) );
        bestllh(i) = matches(1);
    end
end

save bestllh.out bestllh -ASCII -tabs
save llh.out llh -ASCII -tabs

% save parSmean.out par.Smean -ASCII -tabs
% for i = 1:nclusopt
%     save strcat('parSprec', num2str(i), '.out') par.Smean{i} -ASCII -tabs
% end
% save parFmean.out par.Fmean -ASCII -tabs


expect = getExpectations( par );
figure; plot(par.Fmean'); title('Inferred Fmeans');
%figure; plot(par.mumean); title('Inferred mu means');

postD = par.Fmean(:,1);
for i = 2:timepoints
    postD = [ postD par.Smean*postD(:,i-1)+par.Sconstmean ];
end
figure; plot(postD'); title('Inferred Dynamics using S and par.Fmean(:,1)');

max(max(par.Sprec{2}))

if( size(preD) == size(par.Fmean) )
    [ ab ind ] = sort(par.Fmean(:,1));
    sortedF = par.Fmean(ind,:);
    sortedPostD = postD(ind,:);
    [ ab ind ] = sort(preD(:,1));
    sortedPreD = preD(ind,:);
    sortedF - sortedPreD
    sortedPostD - sortedPreD

    [ ab ind ] = sort(diag(par.Smean));
    sortedS = par.Smean(ind,ind);
    [ ab ind ] = sort(diag(Smat));
    sortedSmat = Smat(ind,ind);
    sortedS - sortedSmat
end

%% test
% t = 0;
% 
% simdata;
% 
% for i = 1:10
%     % par = updateDELTA( t, D, par, priors, info );
%     for t = 0:2
%         par = updateMU( t, D, par, priors, info );
%     end
%     % par = updateGAM( t, D, par, priors, info );
%     % par = updateL( t, D, par, priors, info );
% 
%     % for t = info.tuniq
%     %     par = updateF( t, D, par, priors, info );
%     % end
% 
%     % par = updateK( t, D, par, priors, info );
%     % par = updateS( t+1, D, par, priors, info );
% 
%     par = stepClustering(par, info);
% end
