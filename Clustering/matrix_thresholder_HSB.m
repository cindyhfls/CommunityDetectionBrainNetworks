function [kdenth,rth] = matrix_thresholder_HSB(rmat0,th,params,writepajek)    

% modified the script such that the writing of pajek file happens within
% this function such that each iteration just need to add more edges
% instead of having to threshold each time

if ~exist('writepajek','var')||isempty(writepajek)
    writepajek = 1; % default is to write pajek file
end

type = params.type;
writepath=params.writepathbase;

Nroi=length(rmat0); 
NPE=Nroi*(Nroi-1)/2;

% binarize matrix
if params.binary, rmat0=single(rmat0>0); end

switch type
    case 'mst'
         [MST] =backbone_wu(rmat0); %N.B.: I modified the BCT backbone function so that it only gets the backbone itself now
end

[kdenth,rth] = deal(NaN(size(th)));
rmat = triu(rmat0);
switch type
    case 'r'
        [x,y,z] = find(rmat);       % get edges and values
        v = [x,y,z]; v= sortrows(v,3, 'descend' ); 
        ind = arrayfun(@(j)find(v(:,3)>=j,1,'last'), th);
        rth = th;
        kdenth = ind/NPE;
    case 'kden'  
        [x,y,z] = find(rmat);       % get edges and values
        v = [x,y,z]; v= sortrows(v,3, 'descend' );
        ind=ceil(th*NPE);
        rth =v(ind,3);
        kdenth = ind/NPE;
    case 'mst'
        assert((Nroi-1)/NPE<=min(th),'make sure the lowest density is at least (Nroi-1)/NPE'); % 
        notintree = rmat.*(MST==0);
        intree = rmat.*(MST>0);
        [x,y,z] = find(intree); % get edges and values
        v1 = [x,y,z]; v1= sortrows(v1,3, 'descend' );
        [x,y,z] = find(notintree); % get edges and values
        v2= [x,y,z]; v2= sortrows(v2,3, 'descend' );
        v = [v1;v2];
        ind=ceil(th*NPE);
        rth = v(ind,3); % N.B. This is not well defined for the lowest thresholds
        kdenth = ind/NPE;
end

%% Write Pajek file
if writepajek
    for j = 1:length(ind)
        Nedges = ind(j);
        nodenum = 1:Nroi;
        c=clock;
        fprintf(['\t%2.0f:%2.0f:%2.0f: mat2pajek: writing .net file,',...
            'with %d vertices and %d edges\n',c(4),c(5),c(6),Nroi,Nedges]);
        pajekfname = ['paj_col',num2str(j),'.net'];
        outputname=fullfile(writepath,pajekfname);
        fid=fopen(outputname,'wt');
        fprintf(fid,'*Vertices %d\n',Nroi);
        fprintf(fid,'%d "%d"\n',[nodenum; nodenum]);
        fprintf(fid,'*Edges %d\n',Nedges);
        fprintf(fid,'%d %d %f\n',v(1:ind(j),:)');
        fclose(fid);
    end
end
