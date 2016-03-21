% ------------------------------------------------------------------------- 
%                               Constantin's Workflow 
% 
% 
% 
% 
% 
%                                          (c)Constantin Heck, 28-Apr-2015 
%                                                http://mic.uni-luebeck.de
% ------------------------------------------------------------------------- 

clc;
close all;
clear


%what to show
setupPDE           = 1;
makePerfusion      = 0;
makeLocalPerfusion = 0;
saveLargeData      = 0;
saveSmallData      = 0;
dynamicFlow        = 0;

%list of variables to save
varSave = {'Cmat','timeline','aifval','k','n','Hd','qmat','prm','phimat','Fmat','pmat','h','m'};



%configuration
% aif = 'delta';
aif = 'gamma';
% aif = 'parker';
% aif = 'test';
% aif = 'delta2';


%setup spatial-temporal domain
m        = [64,64,1];       %matrix size
omega    = [0,3,0,3,0,3];   %domain in mm
Tmax     = 90;              %max time in seconds
dt       = 0.001;           %time resolution of simulation
tSamp    = .1;              %sampling rate (in seconds)



%setup the other variables
timeline = (0:dt:Tmax); %timeline in seconds
k        = numel(timeline);
h        = (omega(2:2:end)-omega(1:2:end))./m; %grid-size
hd       = prod(h); %voxel-volume, mm^3
Hd       = prod((omega(2:2:end)-omega(1:2:end))); %volume of the whole domain
n        = prod(m);
mk       = [m,k];


%setup basic parameters for the flow
cso  = [1,1,1];
csi  = [64,64,1];
F    = 50/100/60*Hd; %absolute in-/outflow in mm^3/s
K    = 5e-6;         %Permeability: Units in mm^2
mu   = 5*1e-3;       %viscosity is approx 5 cP = 5e-3 Pa*s

%setup basic parameteres for the transport
phimat = .05*ones(m);

%setup prm
prm.cso      = cso;
prm.csi      = csi;
prm.timeline = timeline;
prm.dt       = dt;
prm.scaling  = 1;
prm.mu       = mu;
prm.h        = h;
prm.dim      = m;
prm.physdim  = (omega(2:2:end)-omega(1:2:end-1));


%% setup PDE Model
if setupPDE
    
    
    %setup the flow
    %in-/outflow
    Fmat = zeros(m); 
    for i = 1:size(cso,1)
        Fmat(cso(i,1),cso(i,2),cso(i,3)) = F/size(cso,1);
    end
    for i = 1:size(csi,1)
        Fmat(csi(i,1),csi(i,2),csi(i,3)) = -F/size(csi,1);
    end

    
    %permeability
    Kmat          = K*ones([m,3]); %Kmat is a (diagonal) tensor

    %get the ultimate flow
    [pmat,qmat] = perfusion1c.syntTPFA(Kmat,Fmat,prm); %qmat is absolute flow, i.e. qmat = -K/mu * \nabla p * dA
    pmat        = pmat - min(pmat(:));
    

    %setup the transport

    %setup AIF
    switch aif
        case 'delta'
            aifval    = zeros(size(timeline));
            aifval(1) = phimat(1);

        case 'gamma'
            aifval = perfusion1c.getGammaAIF(timeline/60); %the AIF is assumed to be in mMol/l.
            aifval = aifval*1e-6;                          %Conversion to mMol/mm^3 by multiplication with 1e-6
            
        case 'parker'
            aifval = perfusion1c.getParkerAIF(timeline/60); %the AIF is assumed to be in mMol/l.
            aifval = aifval*1e-6;                          %Conversion to mMol/mm^3 by multiplication with 1e-6
            
        case 'test'
            aifval = perfusion1c.getGammaAIF(timeline/60,'a',2,'b',2); %the AIF is assumed to be in mMol/l.
            aifval = aifval*1e-6;
            
        case 'delta2'
            aifval        = zeros(size(timeline));
            aifval(1:500) = 1;
            aifval        = aifval*1e-6*phimat(1);
    end
    


    %run forward problem    
    disp('Starting simulation...')
    Cmat = perfusion1c.syntforwprobpde(phimat,qmat,Fmat,aifval,prm);
    disp('....done!');

end


%% make perfusion

if makePerfusion    
    [perfmat,lenmat] = perfusion1c.flux2perfCMS(qmat,m,h);
    
    %also save perfusion
    varSave{end+1}   = 'perfmat';
end


%% make local perfusion

if makeLocalPerfusion    
    lperfmat = perfusion1c.localPerfusion(qmat,Fmat,h);
    
    %also save perfusion
    varSave{end+1}   = 'lperfmat';
end



%% save data
if saveLargeData
    
    %save large data
    save('largeDataSet.mat',varSave{:},'-v7.3');
    
end


if saveSmallData
    
    %downsample the data to approx tSamp time resolution
    step     = round(tSamp/dt);
    Cmat     = Cmat(:,:,:,1:step:end);
    timeline = timeline(1:step:end);
    aifval   = aifval(1:step:end);
    k        = numel(timeline);
    
    %save small data    
    save('smallDataSet.mat',varSave{:});
    
end






 %% show the flow

if dynamicFlow

    %display the two indicator flows
    cmax = .1*max(Cmat(:));

    figure(1);clf;
    colormap gray;
    set(1,'name','Idicator-Flow')
    
    for i = 1:k

        Ci  = squeeze(Cmat(:,:,1,i));

        imagesc(Ci);
        axis image
        caxis([0,cmax]);
        ti = sprintf('PDE-Indicator, %1.0f%% completed',i/k*100);
        title(ti);
        
        pause(.1);

        drawnow;

    end


end


