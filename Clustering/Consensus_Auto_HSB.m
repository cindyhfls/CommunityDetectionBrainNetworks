function [stats]=Consensus_Auto_HSB(stats0,cmap,IntMin)
%
% This function takes as input a matrix of cluster assignments (roiXkden),
% and outputs Consensus sortings of the data that explain the primary
% groupings within the data. cmap is the colormap for the Clusters data.
% Goals:
%   1- Minimize parameters and models required to describe data.
%   2- Discard noisy solutions - solutions exhibited at 1 or 2 kden that
%       are not apparent in neighboring solutions.
%   3- Increase the SNR of models that describe data over a range of kden.
% Algorithm:
% 0 Calculates the NMI between each successive column of input. 
% 1 A threshold is iteratively chosen, and neighboring Dkden NMI points are
%   grouped together if they are all above the threshold and if the group 
%   size is above a minimum(parameter: IntMin). Thresholds are chosen
%   following the unique values of the NMI for neighboring kden. This way,
%   all comparisons as th is adjusted are 1-v-1 or 1-v-2.
% 2 A metric of Consensus quality is calculated as the NMI or the Jaccard
%   between all kden sortings within that group. 
% 3 The NMI threshold is then adjusted, groups are found, metrics are
%   calculated. 
% 4 Groups containing the same Dken points as the previous threshold have
%   their metrics compared. 
%       A-  If multiple groups over overlapping kden as a previous th
%           have a metric greater than the single larger group,
%           then they are chosen. 
%       B-  If they are equal,the larger group is preferred as it is a 
%           simpler model for the data set.
%       C-  If the new (higher) NMI threshold provides 1 group with a
%           smaller range of kden, it is chosen over the larger group if
%           the group metric is higher.
% Output- stats structure containing fields:
%       clusters-   input data
%       nmi-        neighboring kden NMI
%       epochs-     final ranges of consensus with group NMI and Jaccard

%% Parameters
Clusters=stats0.SortClus;
[Nroi,Nkden]=size(Clusters);
stats.clusters=double(Clusters);
clear Clusters
if ~exist('IntMin','var'), IntMin=10;end
type=3; % For NMI. to use Jaccard, set type=4.
stats.nmi=zeros(Nkden,1);


%% Calc NMI btwn neighboring kden
for j=1:(Nkden-1)
    foo=nmi_HSB(stats.clusters(:,j+1),stats.clusters(:,j));
    stats.nmi(j)=foo.NMI;
end
th=unique(stats.nmi);
if any(th==0), th(th==0)=[];end
Nth=length(th);
epochs=cell(Nth,1);         % th by th list

%% Use basic thresholding
for k=1:Nth
    disp(['Calc consensus metrics at NMI value = ',num2str(th(k))])
    Ath=find(stats.nmi>=th(k));        % neighboring kden with NMI >= th
    if ~any(Ath)
        continue
    else
        Dath=diff(Ath);             % spacing btwn must be == 1 to combine
        Kf=find(Dath~=1);           % Define segment endpoints
        Ki=cat(1,1,Kf+1);           % set initial endpoints for th
        Kf=cat(1,Kf,length(Dath)+1);
        epochs{k}=[Ki,Kf];          % Cell for th, double array for groups
        keep=find((Kf-Ki+1)>=IntMin);
        epochs{k}=epochs{k}(keep,:);
        NgroupsK=size(epochs{k},1);
        for j=1:NgroupsK
            % Translate epochs back to kden
            epochs{k}(j,1)=Ath(epochs{k}(j,1));
            epochs{k}(j,2)=Ath(epochs{k}(j,2))+1;
            foo=stats.clusters(:,(epochs{k}(j,1)):epochs{k}(j,2));
            foob=Consensus_Metrics_HSB(foo);
            epochs{k}(j,3)=foob.nmi;
            epochs{k}(j,4)=foob.jaccard;
            epochs{k}(j,5)=th(k);
        end
    end  
end


%% Compare epoch metrics against those of previous th
stats.epochs=epochs{1};     % to be master output list
for k=2:Nth
    foo=stats.epochs;
    n=0;                  % If anything splits, this will fix indexing
    NgroupsKm1=size(foo,1);
    NgroupsK=size(epochs{k},1);
    for j=1:NgroupsKm1              % Loop over current master consensus
        rangeKm1=foo(j,1):foo(j,2); % epoch range for current master consensus
        keepK=zeros(NgroupsK,1);
        for l=1:NgroupsK
            rangeK=epochs{k}(l,1):epochs{k}(l,2);       % epoch range for k
            keepK(l)=sum(ismember(rangeKm1,rangeK))>0;
        end
        NkeepK=sum(keepK);
        Dken=find(keepK==1);
        
        % Update the master list
        if NkeepK==0            % No Overlap
            n=n+1;
            stats.epochs(n,:)=foo(j,:);
        elseif NkeepK==1        % New single epoch.
            if foo(j,type)>epochs{k}(Dken,type)
                n=n+1;
                stats.epochs(n,:)=foo(j,:);
            else
                n=n+1;
                stats.epochs(n,:)=epochs{k}(Dken,:);
            end
        else                    % New 2 possible epochs
            Kmetric=0.5*(epochs{k}(Dken(1),type)+...
                epochs{k}(Dken(2),type));
            if foo(j,type)>Kmetric
                n=n+1;
                stats.epochs(n,:)=foo(j,:);
            else
                n=n+1;
                stats.epochs(n,:)=epochs{k}(Dken(1),:);
                n=n+1;
                stats.epochs(n,:)=epochs{k}(Dken(2),:);
            end
        end
    end
end


%% Set consensus of optimized metric epochs
Nc=size(stats.epochs,1);
for j=1:Nc
    foo=stats.clusters(:,(stats.epochs(j,1)):stats.epochs(j,2));
    foo(foo==0) = NaN;% update 2023.09.28 JCT to ignore isolated assignments and use the solutions with an assignment
    tmp = mode(foo');tmp(isnan(tmp)) = 0;
    stats.Consensus(:,j) = tmp;
%     stats.Consensus(:,j)=mode(foo');
end


%% Visualize
Nc=size(stats.Consensus,2);
figure('Color','w','Position',[5,70,1730,1020]);
subplot(2,Nc+1,[1,Nc+2]);
plot(stats.nmi(1:(end-1)),'k');axis([0,Nkden+1,th(1)-0.1,1])
ylabel('NMI of neighboring kden');xlabel('Dkden')
title([{['NMI min = ',num2str(th(1))]};...
    {['Group size min = ',num2str(IntMin),'.']}])


SortRowVec=[1:Nkden,1:Nkden,1:Nkden];
SortRowVecC=[1:Nc,1:Nc,1:Nc];


for j=1:Nc
    mu=round(mean(stats.epochs(j,1:2)));
    subplot(2,Nc+1,j+1);        % full matrix, resorted
    
    foo=stats.clusters;
    foo(foo==0)=NaN;
    SRV=circshift(SortRowVec,[(mu-1),0]);
    for k=1:(3*Nkden)
        foo=sortrows(foo,SRV(k));
    end
    image(sortrows(foo,mu))
    set(gca,'YTick',[]);
    hold on;plot([ones(1,2).*stats.epochs(j,1)],[1,Nroi],'k')
    plot([ones(1,2).*stats.epochs(j,2)],[1,Nroi],'k')
    colormap(cmap)
    freezeColors_HSB
    title([{['Sorted at']};{['kden=',...
        num2str(stats0.kdenth(mu)*100,'%0.2f'),'%']}]);
    xlabel([{['InfoMap']},{['model number']}])
        
    subplot(2,Nc+1,Nc+j+2);     % Consensus matrix, resorted
    foob=stats.Consensus;
    foob(foob==0)=NaN;
    SRVc=circshift(SortRowVecC,[(j-1),0]);
    for k=1:(3*Nc)
        foob=sortrows(foob,SRVc(k));
    end
    image(sortrows(foob,j));
    set(gca,'YTick',[]);
    colormap(cmap)
    freezeColors_HSB
    title([{['kden:']};...
        {[num2str(stats0.kdenth(stats.epochs(j,1))*100,'%0.2f'),'-',...
        num2str(stats0.kdenth(stats.epochs(j,2))*100,'%0.2f'),'%']};...
        {['\mu_N_M_I = ',num2str(stats.epochs(j,type),'%0.3f')]}]);
    xlabel([{['Consensus']},{['model number']}])

    
    subplot(2,Nc+1,[1,Nc+2]);
    vals=stats.epochs(j,1):stats.epochs(j,2);
    hold on;plot(vals,ones(1,length(vals)).*stats.epochs(j,5),...
        'b','LineWidth',5);
end

stats.epochs=table([stats.epochs(:,1)],[stats.epochs(:,2)],...
    [stats.epochs(:,3)],[stats.epochs(:,4)],[stats.epochs(:,5)],...
    'VariableNames',{'kden_i','kden_f','nmi','jaccard','nmi_th'});

