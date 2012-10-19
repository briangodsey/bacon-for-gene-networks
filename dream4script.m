addpath("./src")
addpath("./src/octavemath")

warning("off", "Octave:matlab-incompatible");
warning("off", "all");

%% %%%%%%%%%%%%%%%%%%%%
%% data loading parameters
datdir = "./data/DREAM4/"

numcores = 2
numstarts = 10

%% %%%%%%%%%
%% options for 10-gene networks
numgenes = "10";
ncluses = num2cell([5:10]);
usetimeseries = { [1:5],1,2,3,4,5 };


%% %% %%%%%%%%%%%%
%% %% options for 100-gene networks
%% numgenes = "100";
%% ncluses = num2cell([50:10:100]);
%% usetimeseries = { [1:10],1,2,3,4,5,6,7,8,9,10 };



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
datasets = num2cell([1:5]);
%%datasets = num2cell([1:1]);


%% datasetresults = {};

%%datasetresults = dream4runfunction(numgenes,dsn,datdir);

ndat = size(datasets,2);
ncl = size(ncluses,2);
nts = size(usetimeseries,2);

%% numgeneslist = repmat({numgenes},ndat,ncl);
%% dsnlist = repmat(datasets',1,ncl);
%% datdirlist = repmat({datdir},ndat,ncl);
%% numstartlist = repmat({numstarts},ndat,ncl);
%% ncluslist = repmat(ncluses,ndat,1);

numgeneslist = cell(ndat,nts,ncl);
numgeneslist(:,:,:) = numgenes;

dsnlist = cell(ndat,nts,ncl);
for( i = 1:ndat )
  dsnlist(i,:,:) = datasets(1,i);
end

datdirlist = cell(ndat,nts,ncl);
datdirlist(:,:,:) = datdir;

numstartlist = cell(ndat,nts,ncl);
numstartlist(:,:,:) = numstarts;

ncluslist = cell(ndat,nts,ncl);
for( i = 1:ncl )
  ncluslist(:,:,i) = ncluses(1,i);
end

usetimeserieslist = cell(ndat,nts,ncl);
for( i = 1:nts )
  usetimeserieslist(:,i,:) = usetimeseries(1,i);
end


%% for dsn = 1:5
%%   datasetresults = dream4runfunction(numgenes,dsn,datdir);
%% end

%% ac = {eye(5),eye(7),eye(12)};
%% ab  = parcellfun(3,@size,ac,"UniformOutput",false);

%% %% for debugging, etc, single run
numgenes = numgeneslist{1,1,1};
dsn = dsnlist{1,1,1};
datdir = datdirlist{1,1,1};
nclus = ncluslist{1,1,1};
usets = usetimeserieslist{1,2,1};
shuffle = 0;

%% %%%%%%%%%%%%%%%%%%%%%%%%
%% this runs the model fit
datasetresults = parcellfun(numcores,@dream4runfunction, ...
			    numgeneslist, datdirlist, dsnlist, 
			    usetimeserieslist, ...
			    numstartlist, ncluslist, ...
			    "UniformOutput",false);


%%save finalresults.mat
save(strcat("finalresults_",numgenes,".mat"));


