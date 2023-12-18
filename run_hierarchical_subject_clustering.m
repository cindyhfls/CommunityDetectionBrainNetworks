% replicating Akiki & Abdallah, 2019
% it uses genlouvain but we could switch that out for infomap or other
% models later

listname='washu120' %'washu120'
parcel_name ='Gordon'
ptseriesdir =[ '/data/wheelock/data1/datasets/WashU120/ptseries/']%eLABE/ptseries/Y3/'; % change this
niter = 100;

%'eLABE_Y0_healthy_term_N261'%'eLABE_Y2_N113_atleast600frames'% change this

cohortfile = dir(['/data/wheelock/data1/people/Cindy/BCP/ParcelCreationGradientBoundaryMap/cohortfiles/cohortfiles_',listname,'.txt']);
[subjects, cifti_files, surfdirL, surfdirR] = textread(fullfile(cohortfile.folder,cohortfile.name),'%s%s%s%s');
tmaskfile = dir(['/data/wheelock/data1/people/Cindy/BCP/ParcelCreationGradientBoundaryMap/tmasklist/tmasklist_',listname,'.txt']);
[tmasksubjects, tmaskfiles]=textread(fullfile(tmaskfile.folder,tmaskfile.name),'%s%s');
if ~isequal(tmasksubjects,subjects)
    error('tmasklist subjects do not match cohortfile subjects');
end
all_ci = cell(length(subjects),1);

stats.params.outputdir = '/data/wheelock/data1/people/Cindy/BrBx-HSB_infomap_cleanup/Results/Hierarchical_subject_clustering_Akiki_Abdallah';

%% Obtain community organization for each subject
parfor i =1:length(subjects)
        disp(subjects{i})
        tic
        %% state paths for subject
        [~,subpath] = fileparts(cifti_files{i});
        ptseries_path = fullfile(ptseriesdir,parcel_name,strrep(subpath,'dtseries','ptseries.nii'));
        %% get censored ptseries
        tmask = load(tmaskfiles{i}); % saved as txt
        timeseries = ft_read_cifti_mod(ptseries_path);
        %% get community organization for individual subjects
        all_ci{i,1} = RMT_com(timeseries.data(:,tmask==1)',niter);% 
        toc
end
save(fullfile(stats.params.outputdir,'all_ci_washu120.mat'),'all_ci')
%% (optional) reassign singletons
% for i=1:nindiv
%     tmptmp_ci=all_ci{i,1};
%     newtmp_ci=ci_restoresingleton(tmptmp_ci);
%     all_ci{i,1}=newtmp_ci;
% end
%%
all_ci_combined = horzcat(all_ci{:});
%% Plot the rsults
[Cons] = HierarchicalConsensus_Jeub(all_ci_combined,0.05); % consensus with Jeub et al. 2018 Scientific Reports

% Assign color to consensus
nameoption = 3;% 1: automatic, 3: using template
templatepath  ='Gordon2017_17Networks.dlabel.nii';% 'Tu_eLABE_Y2_22Networks.nii'
parcelpath = '/data/wheelock/data1/parcellations/333parcels/Parcels_LR.dtseries.nii'
[CWro,Cons.SortConsRO] = assign_network_colors(Cons.SortCons,3,templatepath,parcelpath);

parcel_name ='Gordon'%'eLABE_Y2_prelim_072023_0.75'%'Gordon'% params.parcel_name
load(['Parcels_',parcel_name,'.mat'],'Parcels');

warning('off');
Explore_parcel_levels_HSB(Cons.SortCons,CWro.cMap,Parcels,Cons.levels,fullfile(stats.params.outputdir,'consensus'));
%  Explore_parcel_levels_HSB(Cons.modeCons,CWro.cMap,Parcels,Cons.mean_kdenth);
close all;


