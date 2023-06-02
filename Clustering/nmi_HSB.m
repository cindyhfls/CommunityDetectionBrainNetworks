function stats=nmi_HSB(A,B)
%
% This function calculates various information-theoretical metrics for the
% passed variables including: H-entropy, Mab-mutual information,
% NMI-normalized mutual information.
%

%% 
N=numel(A);
A=A(:);
B=B(:);

l=min(min(A),min(B));
A=A-l+1;
B=B-l+1;
u=max(max(A),max(B));

idx=[1:N]';
Ma=sparse(idx,A,1,N,u,N);
Mb=sparse(idx,B,1,N,u,N);


%% Prob Distributions
Pa=mean(Ma,1);          % Individual distributions
Pb=mean(Mb,1);
Pab=nonzeros(Ma'*Mb/N); % Joint distribution


%% Entropies
Ha=-dot(Pa,log2(Pa+eps));
Hb=-dot(Pb,log2(Pb+eps));
Hab=-dot(Pab,log2(Pab+eps));


%% Mutual Information
Mab=Ha+Hb-Hab;
NMI=sqrt((Mab/Ha)*(Mab/Hb));

%% outputs
stats.NMI=NMI;
stats.Mab=Mab;
stats.Ha=Ha;
stats.Hb=Hb;
stats.Pa=Pa;
stats.Pb=Pb;