%% Supporting functions
% FileNameAndLocation=[mfilename('fullpath')];
function mfilebackup(FileNameAndLocation,identifier)
if ~exist('identifier','var')
    identifier = '';
end
timestr=datestr(datetime('now'),'yyyymmdd HH:MM:SS');
[filepath,name] = fileparts(FileNameAndLocation);
if ~exist(fullfile(filepath,'Backup'),'dir')
    mkdir(fullfile(filepath,'Backup'));
end
newbackup=fullfile(filepath,'Backup',sprintf('%s_%s_%s.m',name,timestr,identifier));
currentfile=strcat(fullfile(filepath,[name, '.m']));
copyfile(currentfile,newbackup);
end
