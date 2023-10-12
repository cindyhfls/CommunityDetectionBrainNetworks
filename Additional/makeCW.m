function [CW,GenOrder] = makeCW(params,Nnets)
% manual assignment of network name and colors

switch params.IMap_fn
    case 'Infomap_BCP_220601'
        CW.cMap=zeros(Nnets-1,3);       CW.Nets=cell(Nnets-1,1);n=1;
        CW.cMap(n,:)=[0,1,1];           CW.Nets{n,1}='SMhand';n=n+1;
        CW.cMap(n,:)=[0,0.4,0.4];  CW.Nets{n,1}='VentralAttn';n=n+1;
        CW.cMap(n,:)=[0,1,0];       CW.Nets{n,1}='DorsalAttn';n=n+1;
        CW.cMap(n,:)=[1,0.7,0.9];           CW.Nets{n,1}='Premotor';n=n+1;
        CW.cMap(n,:)=[0,0,1];       CW.Nets{n,1}='VisualLat';n=n+1;
        CW.cMap(n,:)=[1,0,0];       CW.Nets{n,1}='DefaultPost';n=n+1;
        CW.cMap(n,:)=[1.0,0.8,0.6];       CW.Nets{n,1}='VisualMed';n=n+1;
        CW.cMap(n,:)=[1,0.4,0];       CW.Nets{n,1}='SMmouth';n=n+1;
        CW.cMap(n,:)=[1,1,0];       CW.Nets{n,1}='FrontalParietal';n=n+1;
        CW.cMap(n,:)=[0.4,0,0.5];       CW.Nets{n,1}='CinguloOperc';n=n+1;
        CW.cMap(n,:)=[0.1,0.1,0.1];       CW.Nets{n,1}='Salience';n=n+1;
        CW.cMap(n,:)=[1,0.5,0.5];       CW.Nets{n,1}='DefaultAnt';n=n+1;
        CW.cMap(n,:)=[0,0.7,0.7];       CW.Nets{n,1}='VentralAttn2';n=n+1;
        CW.cMap(n,:)=[1,1,0.5];       CW.Nets{n,1}='FrontoParietal2';n=n+1;
%         CW.cMap(n,:) = [0.25,0.25,0.25];   CW.Nets{n,1} = 'None';n = n+1;
        GenOrder = [4,1,3,5,7,12,14,9,8,10,11,2,6,13,14:Nnets-2]; 
    otherwise
        error('cannot find manual color assignment')
end

end
%%
%%
% CW.cMap=zeros(Nnets-1,3);       CW.Nets=cell(Nnets-1,1);n=1;
% CW.cMap(n,:)=[0,0,1];           CW.Nets{n,1}='Vis';n=n+1;
% CW.cMap(n,:)=[1,0,0];           CW.Nets{n,1}='DMN';n=n+1;
% CW.cMap(n,:)=[0,1,1];           CW.Nets{n,1}='Mot';n=n+1;
% CW.cMap(n,:)=[.7,0,.65];         CW.Nets{n,1}='SalN';n=n+1;
% CW.cMap(n,:)=[0.7,0.6,0.2];           CW.Nets{n,1}='SubC';n=n+1;
% 
% CW.cMap(n,:)=[1,1,0];           CW.Nets{n,1}='FPN';n=n+1;
% CW.cMap(n,:)=[1,.5,0];           CW.Nets{n,1}='MotM';n=n+1;
% CW.cMap(n,:)=[0,1,0];         CW.Nets{n,1}='DAN';n=n+1;
% CW.cMap(n,:)=[0,.5,0.15];         CW.Nets{n,1}='MidT';n=n+1;
% CW.cMap(n,:)=[0,0.8,0.7];       CW.Nets{n,1}='VAN';n=n+1;
% 
% 
% CW.cMap(n,:)=[1,1,1];     CW.Nets{n,1}='AUD';n=n+1;
% CW.cMap(n,:)=[0,0.15,0.3];       CW.Nets{n,1}='Mem';n=n+1;
% CW.cMap(n,:)=[0.6,0,0];         CW.Nets{n,1}='DMN2';n=n+1;
% 

% CW.cMap(n,:)=[1,0,1];           CW.Nets{n,1}='COax';n=n+1;
% CW.cMap(n,:)=[0.4,0,0.4];       CW.Nets{n,1}='pCOx';n=n+1;
% 
% CW.cMap(n,:)=[0.4,0.4,0];       CW.Nets{n,1}='FOC';n=n+1;
% CW.cMap(n,:)=[1,0.7,0.7];       CW.Nets{n,1}='aDMN';n=n+1;
% CW.cMap(n,:)=[0.6,0.9,0.4];     CW.Nets{n,1}='daPFC';n=n+1;
% CW.cMap(n,:)=[0.5,0.5,0.5];     CW.Nets{n,1}='Salx';

% CW.cMap(n,:)=[0.2,0.6,0.6];     CW.Nets{n,1}='SMN3';n=n+1;
% CW.cMap(n,:)=[0.4,0.4,0];       CW.Nets{n,1}='FOC';n=n+1;
% CW.cMap(n,:)=[1,0.7,0.7];       CW.Nets{n,1}='aDMN2';n=n+1;
% CW.cMap(n,:)=[0.6,0.9,0.4];     CW.Nets{n,1}='daPFC';n=n+1;
% CW.cMap(n,:)=[0.1,0.1,0.1];     CW.Nets{n,1}='Sal';

%}