function[ datasetresults ] = dream4runfunction(numgenes, datdir, dsn, usets, ...
					       numstarts, nclus)

  disp( strcat("\n\nStarting run: ",numgenes,"_",num2str(dsn), "_", ...
	num2str(usets),"_",num2str(nclus),"\n\n") )

  datasetnum = num2str(dsn);
  %%datasetnum = "1";


  %% %%%%%%%%%%%%%%%%%%%%%
  %% load the gold standard file

  %% goldstdname = strcat("GoldStandardMat_",numgenes,"_", ...
  %% 		       datasetnum,".tsv");
  %% 
  %% disp([goldstdname,"\n"]);
  %% 
  %% goldstd = dlmread(goldstdname)';


  %% %%%%%%%%%%%%%%%%%%%%%%
  %% read in the data file
  filename = strcat(datdir,"DREAM4_InSilico_Size",numgenes, ...
		    "/insilico_size",numgenes,"_", datasetnum, ...
		    "/insilico_size",numgenes,"_",datasetnum,"_timeseries.tsv");
  alldatfile = dlmread(filename);


  %% %%%%%%%%%%%%%%%%%%%%%%%%
  %% pull data set from the loaded matrix
  if( str2num(numgenes) == 10 )
    indALL = {2:(2+20), 23:(23+20), 44:(44+20), ...
	      65:(65+20), 86:(86+20) };
    ind = indALL(usets);
    %%clusternumbers = [ 8 9 10 ];
    %%numstarts = 10;
  elseif( str2num(numgenes) == 100 )
    indALL = {2:(2+20), 23:(23+20), 44:(44+20), ...
	      65:(65+20), 86:(86+20), 107:(107+20), ...
	      128:(128+20), 149:(149+20), 170:(170+20), ...
	      191:(191+20) };
    ind = indALL(usets);
    %%clusternumbers = [ 50 60 70 80 90 100 ];
    %%numstarts = 10;
  end
    
  D = {};
  for i = 1:length(ind)
    D{i} = alldatfile(ind{i},2:end);
  end

  %% %%%%%%%%%%%%
  %% model fit parameters
  %%exper = repmat( 1, 1, 1 );
  names = { "DREAM4_1", "DREAM4_2" "DREAM4_3", ...
	   "DREAM4_4","DREAM4_5" };
  timeORIG = {};
  for i = 1:length(ind)
    timeORIG{i} = alldatfile(ind{i},1)';
  end


  %% %% run parameters
  %% doclustering = 1;
  %% if ~doclustering
  %%   clusternumbers = size(D{1},2);
  %% end



  %% scale-loc normalize the data
  Dorig = D;
  for j = 1:size(D,2)
    for i = 1:size(D{j},2)
      coldat = D{j}(:,i);
      datsd = std(coldat);
      D{j}(:,i) = (coldat-mean(coldat)) / std(coldat);
    end
  end


  infor = setinfo(D,timeORIG,names);

  %% Does the model fitting for various cluster numbers
  %%hyperparvals = 0.0001; %(1:20:200) / ( max(max(D))-min(min(D)) );
  llh = repmat( -1e20, numstarts, 1 ); %%size(clusternumbers,2) );
  %%for i = 1:size(clusternumbers,2)
  %% nclus = clusternumbers(i);
  %% nclus
  
  infor.nclus = nclus;
  priors = setpriors( D, infor );
  
  for j = 1:numstarts
    tic
    [ allpars{j} allpriors{j} llh(j) iter ] = vbmodelfit( D, priors, infor );
    toc
    llhshow = llh(j)
    save(strcat("runresults_",numgenes,"_",datasetnum,".mat"));
  end
  %%end

  
  datasetresults.allpars = allpars;
  datasetresults.allpriors = allpriors;
  datasetresults.llh = llh;
  
  %%save datasetresults.mat
  save(strcat("datasetresults_g",numgenes,"_ds",datasetnum,"_ts", ...
	      num2str(sum(usets)),"_nclus",num2str(nclus),".mat"));
  
end