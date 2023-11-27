function [mods,CodeLength]=infomap_wrapper_HSB(rmat,pajekPath,pajekFn,reps,verbose)
global infomappath
%
% This function generates Pajek files, runs infomap, loads relevant
% modularity assignments to output, and kills Pajek and clu files generated
% by InfoMap.
%% Set some defaults % 2023.09.03
if ~exist('pajekPath','var')||isempty(pajekPath)
    pajekPath = './';
end
if ~exist('pajekFn','var')||isempty(pajekFn)
    pajekFn = 'tmp.net';
end
if ~exist('reps','var')||isempty(reps)
    reps = 1;
end
if ~exist('verbose','var')||isempty(verbose)
    verbose = 0;
end
    
%% Generate Pajek file
here0=pwd;
cd(pajekPath)
outputname=[pajekPath,pajekFn];
if ~exist(outputname,'file')
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
end

%% Obtain seed #
clear randnum;
randnum=ceil(rand*1000000);


%% Run infomap
t=tic;
if verbose
    c=clock;
    fprintf('\t%2.0f:%2.0f:%2.0f: infomap beginning\n',c(4),c(5),c(6));
end

command = [infomappath,' --clu -2 -s',...
    num2str(randnum),' -N',num2str(reps),' ',pajekFn,' ',pajekPath];
[failed, message] = system(command);
if failed
    disp(message)
end

if verbose
    c=clock;
    fprintf('\t%2.0f:%2.0f:%2.0f: infomap finished\n',c(4),c(5),c(6));
    disp(['<<< Infomap took ',num2str(toc(t)),' seconds'])
end

%% Load in modules from clu file
clufile =[pajekPath,pajekFn(1:end-4),'.clu'];
% So parfor doesn't crap out
isclufile = exist(clufile,'file');
while isclufile == 0
    pause(1)
    isclufile = exist(clufile,'file');
end
[mods,CodeLength]=ReadCluFile_HSB(clufile);



%% Clean up files and return
% delete(outputname);
% delete(clufile);
% cd(here0)