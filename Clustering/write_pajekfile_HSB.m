function write_pajekfile_HSB(rmat,pajekPath,pajekFn)
outputname=fullfile(pajekPath,pajekFn);
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