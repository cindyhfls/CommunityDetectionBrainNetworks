function Cons = Cons_stats_HSB(Cons,stats)
%% Stats of consensus models
Cons.stats=Matrix_metrics_HSB(Cons.SortCons,FisherZ2R_HSB(stats.MuMat),...
    Cons.epochs.mean_rth,stats.params.binary);

Ncons=size(Cons.SortCons,2);

%% Plot
figure('position',[100 100 800 800]);
subplot(2,2,1);
plot([1:Ncons],Cons.stats.modularity,'bo-','LineWidth',2);
axis([0,Ncons+1,0,1.1]);
grid on
set(gca,'XTick',[0:Ncons+1])
xlabel('Consensus model');
legend({'Modularity'},'Location','southeast');

subplot(2,2,2);
plot([1:Ncons],Cons.stats.NnBc,'bo-','LineWidth',2);
set(gca,'XTick',[0:Ncons+1])
grid on
xlabel('Consensus model');
legend({'Proportion of nodes in the largest component'},'Location','southeast');

subplot(2,2,3);
plot([1:Ncons],Cons.stats.AvgSil,'bo-','LineWidth',2);
grid on
set(gca,'XTick',[0:Ncons+1])
xlabel('Consensus model');
legend({'Average Silhouette'},'Location','southeast');

subplot(2,2,4);
plot([1:Ncons],Cons.stats.num_communities,'bo-','LineWidth',2);
grid on
set(gca,'XTick',[0:Ncons+1])
xlabel('Consensus model');
legend({'Number of Communities'},'Location','southeast');

