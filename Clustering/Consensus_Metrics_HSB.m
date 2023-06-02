function metrics=Consensus_Metrics_HSB(data)
%
% This function outputs metrics of similarity for the NroiXNkden matrix
% passed in data. Currently supported metrics include: mean NMI and mean
% Jaccard. Basically, these metrics are calculated between each pair of the
% possible kden solutions.

%% Parameters
[Nroi,Nkden]=size(data);
NMIpairs=zeros(Nkden*(Nkden-1)/2,1);
Jacpairs=NMIpairs;
data(data==0)=1000; % switch non-zero values to some number obv not used
n=0;

%% Calculate NMI and Jaccard
for j=1:(Nkden-1)
    foo=squeeze(data(:,j));
    for k=(j+1):Nkden
        n=n+1;
        foob=squeeze(data(:,k));
        temp=nmi_HSB(foo,foob);
        NMIpairs(n)=temp.NMI;
    end
end
metrics.nmi=mean(NMIpairs);
metrics.jaccard=mean(1-pdist(data','jaccard'));