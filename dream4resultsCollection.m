%% addpath('~/fat/research/boku/src/matlab')
%% addpath('~/fat/research/geneInteraction/proj/src')
addpath('./src')

warning ("off", "Octave:matlab-incompatible");

%% %%%%%%%%%%%%%%%%%%%%%%
%% load some saved results


%% %% %%%%%%%%
%% %% old format for results:
%% resultsformat = "oldresults";
%% %% load "~/fat/research/geneInteraction/proj/datasetresults.mat"
%% %% load datasetresults.mat %% contains dsn 1:4
%% %% load finalresults.mat %% don't use
%% load "results10gene_run1/runresults.mat"  %% contains dsn 1:4 in datasetresults{} and dsn 5 in allpars, etc
%% datasetresults{5}.allpars = allpars;
%% datasetresults{5}.allpriors = allpriors;
%% datasetresults{5}.llh = llh;
%% numdatasets = length(datasetresults);


%% %% %%%%%%%%
%% %% new format for results
%% resultsformat = "newresults";
%% %%load "results10gene_run2faster/finalresults_10.mat";
%% %%load "results10gene_run3_10starts/finalresults_10.mat";
%% load "results10gene_run4_10starts_betterllh/finalresults_10.mat";
%% numdatasets = size(datasetresults,1);

%% %%%%%%%%
%% new format for results
resultsformat = "new2results";
%%load "results10gene_run6_5starts_5to10clus_separated/finalresults_10.mat";
%%load "results10gene_run7_5starts_5to10clus_sep_longerrun/finalresults_10.mat";
load "~/fat/research/geneInteraction/runresults/results10gene_FINAL/dream4results_small.mat";
numdatasets = size(datasetresults,1);



numgenes = "10";
%%numgenes = "100";


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
allresults = {};
for dsn = 1:numdatasets

  %%dsn = 5;

  goldstdname = strcat("goldstdmats/GoldStandardMat_",numgenes,"_", ...
		       num2str(dsn),".tsv");
  disp([goldstdname,"\n"]);
  goldstd = dlmread(goldstdname)';


  if strcmp(resultsformat,"oldresults")
    numclusruns = size(datasetresults{dsn}.allpars,2);
    numstarts = size(datasetresults{dsn}.allpars,1);
  elseif strcmp(resultsformat,"newresults")
    numclusruns = size(ncluses,2);
    %%numstarts = size(datasetresults{dsn,1}.allpars,2);
  elseif strcmp(resultsformat,"new2results")
    numclusruns = size(ncluses,2);
    numts = size(usetimeseries,2);
    %%numstarts = size(datasetresults{dsn,1}.allpars,2);
  end

  %% allpars = datasetresults{dsn}.allpars;
  %% allpriors = datasetresults{dsn}.allpriors;
  %% llh = datasetresults{dsn}.llh;
  
  %% par = allpars{ find( llh==max(max(llh)) ) };
  %% priors = allpriors{ find( llh==max(max(llh)) ) };
  %% 
  %% ind = [1,2];
  %% par = allpars{ ind };
  %% priors = allpriors{ ind };
  %% 
  %% [ auroc aupr ] = checkInteractionResults( par, goldstd )
  
  
  aurocmat = [];
  auprmat = [];
  allllhs = [];
  inferredints = {};
  for i = 1:numstarts %%1:size(allpars,1)
    for j = 1:numclusruns %%size(allpars,2)

      if strcmp(resultsformat,"oldresults")
	[ aurocmat(i,j) auprmat(i,j) inferredints{i,j} ] = ...
	    checkInteractionResults( datasetresults{dsn}.allpars{i,j}, goldstd );
	allllhs(i,j) = datasetresults{dsn}.llh(i,j);
      elseif strcmp(resultsformat,"newresults")
	[ aurocmat(i,j) auprmat(i,j) inferredints{i,j} ] = ...
	    checkInteractionResults( datasetresults{dsn,j}.allpars{1,i}, goldstd );
	allllhs(i,j) = datasetresults{dsn,j}.llh(i,1);
      elseif strcmp(resultsformat,"new2results")
	for k = 1:numts
	  [ aurocmat(i,j,k) auprmat(i,j,k) inferredints(i,j,k) ] = ...
	      checkInteractionResults( datasetresults{dsn,k,j}.allpars{1,i}, goldstd );
	  allllhs(i,j,k) = datasetresults{dsn,k,j}.llh(i,1);
	end
      end
    end
  end

  allresults{dsn}.allllhs = allllhs;
  allresults{dsn}.aurocmat = aurocmat;
  allresults{dsn}.auprmat = auprmat;
 
end



resultssummary = [];
resultssummary.auroc = [];
resultssummary.aupr = [];
shortsummary = [];
shortsummary.auroc = [];
shortsummary.aupr = [];
shortsummaryNOCLUS.auroc = [];
shortsummaryNOCLUS.aupr = [];
if strcmp(resultsformat,"oldresults") || strcmp(resultsformat,"newresults")
  for dsn = 1:length(allresults)
    for nc = 1:size(allresults{dsn}.allllhs,2)
      whichbest = find( allresults{dsn}.allllhs(:,nc) ...
		       ==max(allresults{dsn}.allllhs(:,nc)) );
      resultssummary.llh(dsn,nc) = allresults{dsn}.allllhs(whichbest,nc);
      resultssummary.auroc(dsn,nc) = allresults{dsn}.aurocmat(whichbest,nc);
      resultssummary.aupr(dsn,nc) = allresults{dsn}.auprmat(whichbest,nc);
    end
    whichbest = find( resultssummary.llh(dsn,:) ...
		     ==max( resultssummary.llh(dsn,:) ) );
    shortsummary.auroc(dsn) = resultssummary.auroc(dsn,whichbest);
    shortsummary.aupr(dsn) = resultssummary.aupr(dsn,whichbest);
  end
elseif strcmp(resultsformat,"new2results")
  for dsn = 1:length(allresults)
    for nts = 1:numts
      for nc = 1:size(allresults{dsn}.allllhs,2)
	whichbest = find( allresults{dsn}.allllhs(:,nc,nts) ...
			 ==max(allresults{dsn}.allllhs(:,nc,nts)) );
	resultssummary.llh(nts,nc,dsn) = allresults{dsn}.allllhs(whichbest,nc,nts);
	resultssummary.auroc(nts,nc,dsn) = allresults{dsn}.aurocmat(whichbest,nc,nts);
	resultssummary.aupr(nts,nc,dsn) = allresults{dsn}.auprmat(whichbest,nc,nts);
      end
      whichbest = find( resultssummary.llh(nts,:,dsn) ...
		       ==max( resultssummary.llh(nts,:,dsn) ) );
      shortsummary.auroc(nts,dsn) = resultssummary.auroc(nts,whichbest,dsn);
      shortsummary.aupr(nts,dsn) = resultssummary.aupr(nts,whichbest,dsn);

      shortsummaryNOCLUS.auroc(nts,dsn) = resultssummary.auroc(nts,end,dsn);
      shortsummaryNOCLUS.aupr(nts,dsn) = resultssummary.aupr(nts,end,dsn);
    end
  end
end  


save allresults.mat allresults
save resultssummary.mat resultssummary
save shortsummary.mat shortsummary
save shortsummaryNOCLUS.mat shortsummaryNOCLUS





