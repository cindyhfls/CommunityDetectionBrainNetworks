function Plot_InfoMap_Metrics_HSB(stats,clusters)
%
%
[Nroi,Nanals]=size(clusters);

figure('Color','w','Position',[20,110,900,800])

%% kave and ln(N) as f(kden)
tooSparseTh=find(stats.kave>log(Nroi),1);
subplot(2,2,1) 
plot(stats.kdenth,ones(Nanals,1).*log(Nroi),'bd',...
    'Markersize',4,'MarkerFaceColor','b');hold on
plot(stats.kdenth,stats.kave,'rs','Markersize',4,'MarkerFaceColor','r');
x=[min(stats.kdenth),stats.kdenth(tooSparseTh),...
    stats.kdenth(tooSparseTh),min(stats.kdenth)];
y=[0,0,max(stats.kave),max(stats.kave)];
v=[x',y'];
f=[1,2,3,4];
patch('Faces',f,'Vertices',v,'FaceColor',[0.25,0.25,0.25],'FaceAlpha',0.25)
axis([0.5*min(stats.kdenth),max(stats.kdenth),0,max(stats.kave)]);
legend('Ln ( Nroi )','k ave','Location','northwest');legend BOXOFF
xlabel('Edge Density th');

% text(double(min(stats.kdenth)),1,[{['Network may']},{['be too sparse']}])

subplot(2,2,2) 
plot(stats.rth,ones(Nanals,1).*log(Nroi),'bd',...
    'Markersize',4,'MarkerFaceColor','b');hold on
plot(stats.rth,stats.kave,'rs','Markersize',4,'MarkerFaceColor','r');
x=[max(stats.rth),stats.rth(tooSparseTh),...
    stats.rth(tooSparseTh),max(stats.rth)];
y=[0,0,max(stats.kave),max(stats.kave)];
v=[x',y'];
f=[1,2,3,4];
patch('Faces',f,'Vertices',v,'FaceColor',[0.25,0.25,0.25],'FaceAlpha',0.25)
axis([0.9*min(stats.rth),max(stats.rth),0,max(stats.kave)]);
legend('Ln ( Nroi )','k ave','Location','northeast');legend BOXOFF
xlabel('Minimum r-value weights in network model');





%% Modularity and Assortativity
% subplot(3,2,3) 
% plot(stats.kdenth,stats.modularity,'bd',...infinite
%     'Markersize',4,'MarkerFaceColor','b');hold on
% plot(stats.kdenth,stats.A,'rs','Markersize',4,'MarkerFaceColor','r');
% xlabel('Edge Density th');
% axis([0.5*min(stats.kdenth),max(stats.kdenth),0,1]);
% legend('Modularity','Assortativity','Location','northeast');legend BOXOFF
% 
% subplot(3,2,4) 
% plot(stats.rth,stats.modularity,'bd',...
%     'Markersize',4,'MarkerFaceColor','b');hold on
% plot(stats.rth,stats.A,'rs','Markersize',4,'MarkerFaceColor','r');
% xlabel('r th');
% axis([0.5*min(stats.rth),max(stats.rth),0,1]);
% legend('Modularity','Assortativity','Location','northwest');legend BOXOFF


%% Ncomps,and percentage in biggest comp, connectedness
subplot(2,2,3) 
plot(stats.kdenth,stats.Nc,'kp','Markersize',4,'MarkerFaceColor','k');
ax1=gca;ax2=axes;
plot(stats.kdenth,stats.NnBc,'bs','Parent',ax2,...
    'Markersize',4,'MarkerFaceColor','b');hold on
% plot(stats.kdenth,stats.Cdns,'ro','Parent',ax2,...
%     'Markersize',4,'MarkerFaceColor','r');
set(ax1,'box','off')
set(ax2,'Position',get(ax1,'Position'),'box','off',...
    'YAxisLocation','Right','Color','none','XTick',[])
xlabel(ax1,'Edge Density th');
axis(ax1,[0.5*min(stats.kdenth),max(stats.kdenth),0,max(stats.Nc)]);
axis(ax2,[0.5*min(stats.kdenth),max(stats.kdenth),0,1]);
ylabel(ax1,'Number of components');
ylabel(ax2,'%nodes in largest component')
tooDisconnected=min([find(stats.Nc<10,1),find(stats.NnBc>0.9,1)]);
x=[min(stats.kdenth),stats.kdenth(tooDisconnected),...
    stats.kdenth(tooDisconnected),min(stats.kdenth)];
y=[0,0,1,1];
v=[x',y'];
f=[1,2,3,4];
patch('Faces',f,'Vertices',v,'FaceColor',[0.25,0.25,0.25],'FaceAlpha',0.25)

x0=stats.kdenth(tooDisconnected)+0.01*(max(stats.kdenth)-min(stats.kdenth));
y0=0.4;
x=[x0,x0+0.005,x0+0.005,x0];
y=[y0,y0,y0+0.1,y0+0.1];
v=[x',y'];
f=[1,2,3,4];
% patch('Faces',f,'Vertices',v,'FaceColor',[0.25,0.25,0.25],'FaceAlpha',0.25)

text(double(x0+0.001),0.5,[{['kden/r-value range at which']};...
                    {['networks may be too disjointed for ']};...
                    {['accurate interpretation of clustering']}])



subplot(2,2,4)
plot(stats.rth,stats.Nc,'kp','Markersize',4,'MarkerFaceColor','k');
ax1=gca;ax2=axes;
plot(stats.rth,stats.NnBc,'bs','Parent',ax2,...
    'Markersize',4,'MarkerFaceColor','b');hold on
% plot(stats.rth,stats.Cdns,'ro','Parent',ax2,...
%     'Markersize',4,'MarkerFaceColor','r');
set(ax1,'box','off')
set(ax2,'Position',get(ax1,'Position'),'box','off',...
    'YAxisLocation','Right','XTick',[],'Color','none')
xlabel(ax1,'Minimum r-value weights in network model');
axis(ax1,[0.9*min(stats.rth),max(stats.rth),0,max(stats.Nc)]);
axis(ax2,[0.9*min(stats.rth),max(stats.rth),0,1]);
ylabel(ax1,'Number of components');
ylabel(ax2,'%nodes in largest component')
tooDisconnected=min([find(stats.Nc<10,1),find(stats.NnBc>0.9,1)]);
x=[max(stats.rth),stats.rth(tooDisconnected),...
    stats.rth(tooDisconnected),max(stats.rth)];
y=[0,0,1,1];
v=[x',y'];
f=[1,2,3,4];
patch('Faces',f,'Vertices',v,'FaceColor',[0.25,0.25,0.25],'FaceAlpha',0.25)
% text(double(max(stats.rth)),0.2,[{['Network may']},{['be too sparse']}])

