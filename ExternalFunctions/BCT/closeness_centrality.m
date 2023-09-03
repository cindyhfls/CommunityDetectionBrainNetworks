% Created 2022.04.22 Jiaxin Cindy Tu
% N.B. this is custom code not from BCT toolbox
% adapted from weight_conversion.m and distance_wei.m from BCT toolbox
% reference: chapter 5, Fundamentals of Brain Network Analysis
%            Opsahl et al., 2010; Rochat 2009 to calculate closeness
%            centrality for fragmented networks

function Cc=closeness_centrality(W)
    L = W;
    E=find(W);
    L(E)=1./W(E);           % invert weights
    
    n=length(L);
    D=inf(n);
    D(1:n+1:end)=0;                             %distance matrix
    B=zeros(n);                                 %number of edges matrix
    
    for u=1:n
        S=true(1,n);                            %distance permanence (true is temporary)
        L1=L;
        V=u;
        while 1
            S(V)=0;                             %distance u->V is now permanent
            L1(:,V)=0;                          %no in-edges as already shortest
            for v=V
                T=find(L1(v,:));                %neighbours of shortest nodes
                [d,wi]=min([D(u,T);D(u,v)+L1(v,T)]);
                D(u,T)=d;                       %smallest of old/new path lengths
                ind=T(wi==2);                   %indices of lengthened paths
                B(u,ind)=B(u,v)+1;              %increment no. of edges in lengthened paths
            end
            
            minD=min(D(u,S));
            if isempty(minD)||isinf(minD)       %isempty: all nodes reached;
                break,                          %isinf: some nodes cannot be reached
            end;
            
            V=find(D(u,:)==minD);
        end
    end
    Cc = (n-1)./sum(D);
end