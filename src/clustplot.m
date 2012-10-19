function[] = clustplot(X,ctrs,idx)

figure;
plot(X(idx==1,1),X(idx==1,2),'r.','MarkerSize',12)
hold on
plot(X(idx==2,1),X(idx==2,2),'b.','MarkerSize',12)
plot(X(idx==3,1),X(idx==3,2),'g.','MarkerSize',12)
plot(X(idx==4,1),X(idx==4,2),'c.','MarkerSize',12)
%plot(ctrs(:,1),ctrs(:,2),'kx','MarkerSize',12,'LineWidth',2)
plot(ctrs(:,1),ctrs(:,2),'ko','MarkerSize',12,'LineWidth',1)
%legend('Cluster 1','Cluster 2','Centroids','Location','NW')
hold off
end