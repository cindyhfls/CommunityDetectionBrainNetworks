function Cons = Cons_metrics_HSB(Cons,stats)
%% Stats of consensus models
if isfield(Cons,'kdenth')
    levels = Cons.kdenth;
elseif isfield(Cons,'K')
    levels = Cons.K;
elseif isfield(Cons,'gamma')
    levels = Cons.gamma;
end

if ~isfield(Cons,'metrics')||isempty(Cons.metrics)
    Cons.metrics=Matrix_metrics_HSB(Cons.SortCons,FisherZ2R_HSB(stats.MuMat));
end
Ncons=size(Cons.SortCons,2);

%% Plot
figure('position',[100 100 1200 800]);

% subplot(3,2,1);
% plot([1:Ncons],Cons.metrics.num_communities,'bo-','LineWidth',2);
% grid on
% set(gca,'XTick',round(linspace(0,Ncons+1,10)))
% xlabel('Consensus model');
% title('Number of Communities');

subplot(3,2,1);
plot([1:Ncons],Cons.metrics.non_singleton,'b-','LineWidth',2);
grid on
% set(gca,'XTick',round(linspace(0,Ncons+1,10)))
xlabel('Consensus model');
title('Number of Non-singleton Communities');

subplot(3,2,2);
plot([1:Ncons],Cons.metrics.SysSeg,'b.-','LineWidth',2);
grid on
% set(gca,'XTick',round(linspace(0,Ncons+1,10)))
xlabel('Consensus model');
title('System segregation');


subplot(3,2,3);
plot([1:Ncons],Cons.metrics.SilTemporal,'b-','LineWidth',2);
grid on
% set(gca,'XTick',round(linspace(0,Ncons+1,10)))
xlabel('Consensus model');
title('Average Silhouette');
% legend({'Average Silhouette'},'Location','southeast');

subplot(3,2,4);
plot([1:Ncons],Cons.metrics.DB,'b-','LineWidth',2);
% set(gca,'XTick',round(linspace(0,Ncons+1,10)))
grid on
xlabel('Consensus model');
% legend({'DB'},'Location','southeast');
title('Davies-Bouldin');

subplot(3,2,5);
plot([1:Ncons],Cons.metrics.CH,'b-','LineWidth',2);
% set(gca,'XTick',round(linspace(0,Ncons+1,10)))
grid on
xlabel('Consensus model');
% legend({'CH'},'Location','southeast');
title('Calinski-Harabasz');

subplot(3,2,6);hold on;
yyaxis left
plot([1:Ncons],Cons.metrics.Qsigned,'b-','LineWidth',2);
yyaxis right
plot([1:Ncons],Cons.metrics.QNG,'g-','LineWidth',2);
grid on
% set(gca,'XTick',round(linspace(0,Ncons+1,10)))
xlabel('Consensus model');
legend({'Signed','Newman'},'Location','southeast');
title('Modularity')
