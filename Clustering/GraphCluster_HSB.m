function stats=GraphCluster_HSB(rmat,params)

% This function is an adaptation of ngt.graphcluster that does not  
% access the directory structure, but rather is passed it from
% the MATLAB workspace. Similarly, it outputs to the workspace, not  
% the directory structure. Bootstrapping is not supported.
% This program is currently made to be run in the folder containing the
% *.roi file (or at least have it in your path).
% Inputs:
%   rmat            matrix, Nroi-Nroi-Nsubjm, of Pearson correlations
%   params...
%       .binary     binarize at threshold
%       .type       'r' or 'kden'
%       .lo         lower value of rval or kden
%       .step       step size for analyses parameter
%       .hi         high value for rval or kden
%       .xdist      if >0, kill edges where roi sep <= xdist
%       .roi        3D roi coords for xdist flag
%       .roifname   *.roi filename for infomap call 
%               -->(make sure folder with this is in the MATLAB path)
%       .pajekfname  pajek filename for infomap call
%       .writepath  absolute path for output files
%       .fig        to create a figure of basic graph data: 
%               -->all as f(Rth,Kdenth): Nedges(K), edge_density(kden),
%                   average edges per node(kave), modularity(Q), 
%                   assortativity(A), #comps(Nc), %n_in_giant_comp(NnBc)
%
% Outputs:
%   stats.clusters                  cluster solutions
%   stats.params                    params information
%   stats.MuMat                     mean matrix used for infomapping

%   stats.rth                       rth
%   stats.kdenth                    kdenth
%   stats.Nedges                    Nedges     


%% Set up default parameters and initialize stuff
here=pwd;
if ischar(rmat),load(rmat);rmat=corrmat;clear corrmat;end
if ~exist('params','var'), params=struct; end
if ~isstruct(params)            % If params is name of *.prm file, load in
    fname=params;
    fileID=fopen(fname);
    params=struct;
    stuff=textscan(fileID,'%s');                % read the prmfile
    params.prmfile = fileID;                    % name of prm file
    params.roifname=stuff{2,1};                 % name of roi file
    params.lo=str2double(stuff{6,1});           % lo r or kden
    params.step=str2double(stuff{7,1});         % step in r or kden
    params.hi=str2double(stuff{8,1});           % hi r or kden
    params.writepathbase=stuff{9,1};                % write path
    params.type=stuff{13,1};                    % r or kden
    fclose(fileID);
end
if ~isfield(params,'binary'), params.binary=0;end
if ~isfield(params,'type'), params.type='kden';end
if ~isfield(params,'lo'), params.lo=0.01;end % 0.005
if ~isfield(params,'step'), params.step=0.001;end %0.001
if ~isfield(params,'hi'), params.hi=0.10;end % 0.25, 0.1
if ~isfield(params,'xdist'), params.xdist=21;end % was 30, sri used 21
if ~isfield(params,'numworkers'),params.numworkers = 10; end
if ~isfield(params,'writepathbase') % this folder must already exist
    params.writepathbase=[pwd,'/'];
end

if params.xdist                   % exclusion distance
if ~isfield(params,'roi') && isfield(params,'roifname')        % if yes, need coords of ROI centers
    foo=fopen(params.roifname);
    foob=textscan(foo,'%f %f %f %f %s %s %s %s %s','HeaderLines',1);
    params.roi=[foob{1,2},foob{1,3},foob{1,4}];
    fclose(foo);
    clear foo foob
end
end
if ~isfield(params,'fig'), params.fig=1;end

%% assertions
sz = size(rmat);
assert(length(sz)==2,'input has to be only 2-D'); 
assert(sz(1)==sz(2),'input has to be squared'); % input has to be squared
assert(sum(isnan(rmat(:)))==0,'input has to not contain NaN') % make sure NaNs have been removed already they should not contribute to the assignment or edge density

%% Initialize outputs
Nroi = sz(1);
th=params.lo:params.step:params.hi;
Nkden=length(th);
mods=zeros(Nroi,Nkden,'single');       % Clustering labels
codelength = zeros(1,Nkden);   
[rth,kdenth]=deal(zeros(Nkden,1,'single'));    % Modularity % r threshold % kden threshold  % Silhouette index           
if params.repeats_consensus
        stats.avgVersatility = NaN(1,Nkden);
        stats.stdVersatility = NaN(1,Nkden);
        stats.Versatility = NaN(Nroi,Nkden);
end
%% Organize file structure to a temp to kill folder tree
cd(params.writepathbase)
here0=pwd;
foobr=ceil(1e7*rand(1));
GCtempFname=['temp',num2str(foobr)];
mkdir(GCtempFname)
params.writepathbase=[params.writepathbase,GCtempFname,'/']; 
cd(params.writepathbase) 
% writepath=params.writepathbase;

%% Distance exclusion
% To-do: use a generative model for soft distance exclusion
if isfield(params,'dmat')
    dmat=params.dmat>params.xdist;
else
    dmat=pdist2(params.roi,params.roi)>params.xdist;
end
rmat=rmat.*dmat;
clear dmat

%% Remove diagonal correlation
rmat0 = rmat;
for i = 1:Nroi
    rmat0(i,i) = 0;
end

%% Preparation
if params.repeats_consensus
    Versatility = NaN(length(rmat),params.reps);
    [avgVersatility,stdVersatility] = deal(NaN(1,params.reps));
end

if isempty(gcp('nocreate')) && params.numworkers>0
    parpool(params.numworkers);
end
[allpartitions,allcodelengths] = deal(cell(1,Nkden));

%% Threshold and write the file to pajek
[kdenth,rth]  = matrix_thresholder_HSB(rmat0,th,params,1);    

%% Run Infomap on pajek
parfor j=1:Nkden    
    disp(['Model ',num2str(j),' of ',num2str(Nkden)])
     
    % call infomap
    if ~params.repeats_consensus
         [mods(:,j),codelength(j)]=run_infomap_on_pajekfile_HSB(pajekfname,writepath,params.repeats,params.infomappath,params.version);
    else
        partitions = NaN(Nroi,params.repeats);
        cl= NaN(1,params.repeats);
        for k = 1:params.repeats
            [partitions(:,k),cl(k)] =run_infomap_on_pajekfile_HSB(pajekfname,writepath,1,params.infomappath,params.version);
        end
        %         [Versatility(:,j),avgVersatility(j),stdVersatility(j)] = get_nodal_versatility(partitions); % average versatility
        % use the highest quality/lowest codelength partition
        codelength(j) = min(cl);
        mods(:,j) = partitions(:,find(cl==min(cl),1));
        allpartitions{j} = partitions;
        allcodelengths{j} = cl;
        
    end
    delete(fullfile(writepath,pajekfname));
end

cd(here0)
rmdir(GCtempFname,'s')
cd(here);

stats.clusters=mods;            
stats.params=params;            
stats.rth=rth;                  
stats.kdenth=kdenth;            
stats.MuMat=rmat0;              
stats.codelength = codelength; 

if params.repeats_consensus
    stats.allpartitions = allpartitions;
    stats.allcodelengths = allcodelengths;
%     stats.Versatility = Versatility;
%     stats.avgVersatility = avgVersatility;
%     stats.stdVersatiltiy = stdVersatility;
end

clear params;
clear mods;
clear rth;   
clear kdenth;  
clear codelength
clear rmat0;
