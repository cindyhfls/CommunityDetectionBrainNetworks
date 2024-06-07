%% If loading from txt
clear
% params.outputdir = '/data/wheelock/data1/people/Cindy/BCP/Infomap/vertex-wise/WashU120/20230502'

% params.outputdir = '/data/wheelock/data1/people/Cindy/BCP/Infomap/vertex-wise/eLABE_Y2/20230509'
params.outputdir = '/data/wheelock/data1/people/Cindy/BCP/Infomap/vertex-wise/eLABE_Y2_N92/20231102'
stats.kdenth = importdata(fullfile(params.outputdir,'thresholds.txt'));
load(fullfile(params.outputdir,'params.mat'));

if ~isfield(params,'removalidx')
    params.removalidx = false(59412,1);
    if exist(fullfile(params.outputdir,'NaN_idx.mat'),'file')
        load(fullfile(params.outputdir,'NaN_idx.mat'),'idx')
        params.removalidx(idx) = true;
    end
end

if sum(params.removalidx)>0
    stats.clusters = zeros(length(params.removalidx),size(stats.kdenth,2));
    stats.clusters(~params.removalidx,:) = importdata(fullfile(params.outputdir,'rawassn.txt'));
end

%%
minsize = 400;
stats.SortClus   = postprocess_ordinal_multilayer(stats.clusters);  % sort across level to have some consistency across columns
stats.SortClus = remove_singleton(stats.SortClus ,minsize);
stats.SortClus = remove_few_column_clusters(stats.SortClus);

nameoption = 3;% 1: automatic, 3: using template
templatepath  ='Gordon2017_17Networks.dlabel.nii'%'Sylvester2022_communities.dlabel.nii'%'Gordon2017_17Networks.dlabel.nii';%'Tu_eLABE_Y2_22Networks.nii'%'Gordon2017_17Networks.dlabel.nii';
ParcelCommunities =cifti_read(templatepath); % still use the Gordon colors for an unknown parcel?
[CW,GenOrder] = assign_networks_by_template_cifti(stats.SortClus,ParcelCommunities,0.1,'dice')

% The following code prepares the network order and color infomation
CWro.Nets=CW.Nets(GenOrder);
CWro.cMap=CW.cMap(GenOrder,:);
foo=stats.SortClus;foo(foo==0)=NaN;
stats.SortClusRO=stats.SortClus;
for j=1:length(GenOrder),stats.SortClusRO(foo==GenOrder(j))=j;end
load('Parcels_cortex_nomedialwall.mat')
%% (optional) Viewing and Manual edit of specific networks
for iNet =11:length(CWro.Nets)
    Edit_NetworkColors(stats.SortClusRO,CWro,iNet,Parcels);
    minthresh = min(stats.kdenth(any(stats.SortClusRO==iNet)))*100;
    maxthresh = max(stats.kdenth(any(stats.SortClusRO==iNet)))*100;
    text(0.2,1.2,sprintf('%1.2f%%-%1.2f%%',minthresh,maxthresh),'Units','Normalized');
    print(fullfile(params.outputdir,sprintf('Networks%02d',iNet)),'-dtiff','-r300');
%     pause;
    close all;
end

%% Plot a legend for the Networks
N = length(CWro.Nets)
% figure('Units','inches','position',[10 10 2,3]);%[10 10 5 2]
figure('Units','inches','position',[10 10 5,3]);%[10 10 5 2]
h = gscatter(ones(1,N),ones(1,N),CWro.Nets,CWro.cMap,'s',50);
for i = 1:N
    set(h(i),'Color','k','MarkerFaceColor',CWro.cMap(i,:));
end
legend(CWro.Nets,'interpreter','none','FontSize',10,'location','best','Orientation','horizontal','NumColumns',2);
legend('boxoff')
xlim([10,11]);
axis('off')
print(fullfile(params.outputdir,'networks_Legend'),'-dtiff','-r300');
%% Save solution
cd(params.outputdir)
save(['solution_minsize_',num2str(minsize),'.mat'],'CWro','stats');

%% Plot single column
for thresh = [0.25,1.25,1.75,2.75,4,5.75,10,16,19.25]
    icol = find(round(stats.kdenth*100,2)==thresh)
    plot_network_assignment_parcel_key(Parcels, stats.SortClusRO(:,icol),CWro.cMap,CWro.Nets)
    print(fullfile(params.outputdir,sprintf('kden_%1.4f.tif',stats.kdenth(icol))),'-dtiff','-r300');
    close all
end
%%
Explore_parcel_levels_HSB(stats.SortClusRO,CWro.cMap,Parcels,stats.kdenth,fullfile(params.outputdir,'kden'),40);

%% Make video
Make_parcel_kden_Video(stats.SortClusRO,CWro.cMap,Parcels,stats.kdenth,fullfile(params.outputdir,'video'))

%% Save to dlabel file so that the border can be shown 
load('empty_dlabel.mat','data')
% data=ParcelCommunities;
for k = 1:size(CWro.cMap,1)
    data.diminfo{1,2}.maps.table(k+1).key = k;
    data.diminfo{1,2}.maps.table(k+1).rgba=[CWro.cMap(k,:)';1];
    data.diminfo{1,2}.maps.table(k+1).name = CWro.Nets{k};
end
for col = 1:size(stats.SortClusRO,2)
    data.cdata = single(stats.SortClusRO(:,col));
    cifti_write(data,fullfile(params.outputdir,['solution_minsize_',num2str(minsize),'_kden_',num2str(stats.kdenth(col)),'.dlabel.nii']));
end
disp('done')