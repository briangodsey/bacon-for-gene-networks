
%% 10 genes
load ./results10gene/allresults.mat
load ./results10gene/resultssummary.mat
load ./results10gene/shortsummary.mat
load ./results10gene/shortsummaryNOCLUS.mat


%% find out who wins
dfs = shortsummaryNOCLUS.auroc(2:end,:) - shortsummary.auroc(2:end,:);
sgns = sign(dfs)
sum(sum(sgns<0))
sum(sum(sgns>0))
sum(sum(sgns==0))

dfs.auroc = shortsummaryNOCLUS.auroc(2:end,:) - shortsummary.auroc(2:end,:);
dfs.aupr = shortsummaryNOCLUS.aupr(2:end,:) - shortsummary.aupr(2:end,:);


## *restab* is a table in the paper; how many times clustering is better
## than no clustering
levels = [-1 0 1 ];
for aurocstat = 1:length(levels)
  for auprstat = 1:length(levels)
    restab(aurocstat,auprstat) = sum(sum(
					 sign(dfs.auroc)==levels(aurocstat)
					 ...
					 & \
					 sign(dfs.aupr)==levels(auprstat) \
					 ));
  end
end


%% plot(shortsummary.auroc(2:end,:),shortsummary.aupr(2:end,:),"@33", ...
%%      shortsummaryNOCLUS.auroc(2:end,:),shortsummaryNOCLUS.aupr(2:end,:),"@44")

temp = zeros(5,10);
temp(:,[1 3 5 7 9]) = round(100*shortsummary.auroc(2:end,:))/100;
temp(:,[1 3 5 7 9]+1) = round(100*shortsummaryNOCLUS.auroc(2:end,:))/100;


temp = zeros(5,10);
temp(:,[1 3 5 7 9]) = round(100*shortsummary.aupr(2:end,:))/100;
temp(:,[1 3 5 7 9]+1) = round(100*shortsummaryNOCLUS.aupr(2:end,:))/100;


%% look at best number of clusters
bestnumclus = zeros(size(resultssummary.llh,3),size(resultssummary.llh,1));
ncluses = 5:10;
for i = 1:size(resultssummary.llh,3)
  for j = 1:size(resultssummary.llh,1)
    bestnumclus(i,j) = ncluses(find( resultssummary.llh(j,:,i)==max(resultssummary.llh(j,:,i)) ));
  end
end

sum(sum(bestnumclus(:,2:end)==10))