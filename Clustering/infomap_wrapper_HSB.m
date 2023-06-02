function mods=infomap_wrapper_HSB(rmat,pajekPath,pajekFn,reps)
%
% This function generates Pajek files, runs infomap, loads relevant
% modularity assignments to output, and kills Pajek and clu files generated
% by InfoMap.

%% Generate Pajek file
here0=pwd;
cd(pajekPath)
outputname=[pajekPath,pajekFn];
rmat=triu(rmat,1);          % only the upper triangle
nodes = size(rmat,1); % make the input file
[x,y,z] = find(rmat);       % get edges and values
use = [x,y,z];
nodenum = 1:nodes;
c=clock;
fprintf(['\t%2.0f:%2.0f:%2.0f: mat2pajek: writing .net file,',...
    'with %d vertices and %d edges\n',c(4),c(5),c(6),nodes,length(x)]);
fid=fopen(outputname,'wt');
fprintf(fid,'*Vertices %d\n',nodes);
fprintf(fid,'%d "%d"\n',[nodenum; nodenum]);
fprintf(fid,'*Edges %d\n',length(x));
fprintf(fid,'%d %d %f\n',use');
fclose(fid);


%% Obtain seed #
clear randnum;
randnum=ceil(rand*1000000);


%% Run infomap
t=tic;
c=clock;
fprintf('\t%2.0f:%2.0f:%2.0f: infomap beginning\n',c(4),c(5),c(6));
command = ['!/usr/local/infomap/Infomap --clu -2 -s',...
    num2str(randnum),' -N',num2str(reps),' ',pajekFn,' ',pajekPath];
evalc(command);
c=clock;
fprintf('\t%2.0f:%2.0f:%2.0f: infomap finished\n',c(4),c(5),c(6));
disp(['<<< Infomap took ',num2str(toc(t)),' seconds'])


%% Load in modules from clu file
clufile =[pajekPath,pajekFn(1:end-4),'.clu'];
M=ReadCluFile_HSB(clufile);
mods=M(:,2);


%% Clean up files and return
delete(outputname);
delete(clufile);
cd(here0)