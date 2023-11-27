function [mods,CodeLength]=run_infomap_on_pajekfile_HSB(pajekFn,pajekPath,reps,infomappath,version)%(pajekFn,pajekPath,reps,verbose)

%% Obtain seed #
clear randnum;
randnum=ceil(rand*1000000);

%% Run infomap
t=tic;

c=clock;
fprintf('\t%2.0f:%2.0f:%2.0f: infomap beginning\n',c(4),c(5),c(6));


command = [infomappath,' --clu -2 -s',...
    num2str(randnum),' -N',num2str(reps),' ',fullfile(pajekPath,pajekFn),' ',pajekPath];
[failed, message] = system(command)
if failed
    disp(message)
end


c=clock;
fprintf('\t%2.0f:%2.0f:%2.0f: infomap finished\n',c(4),c(5),c(6));
disp(['<<< Infomap took ',num2str(toc(t)),' seconds'])


%% Load in modules from clu file
clufile =fullfile(pajekPath,strrep(pajekFn,'.net','.clu'));

isclufile = exist(clufile,'file');% So parfor doesn't crap out - is this really necessary?
while isclufile == 0
    pause(1)
    isclufile = exist(clufile,'file');
end

[mods,CodeLength]=ReadCluFile_HSB(clufile,version);
delete(clufile);
