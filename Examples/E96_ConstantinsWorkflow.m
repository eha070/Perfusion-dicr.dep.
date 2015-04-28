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

clear;
clc;
close all;


%what to show
dynamicFlow      = 1;
renewCalculation = 1;

%configuration
% aif = 'delta';
aif = 'gamma';


%setup spatial domain
m        = [64,64,1];       %matrix size
omega    = [0,10,0,10,0,1]; %domain in mm
h        = (omega(2:2:end)-omega(1:2:end))./m;
hd       = prod(h); %voxel-volume

%setup timeline
dt       = 0.002;           %in seconds
Tmax     = 90;
timeline = (0:dt:Tmax); %timeline in seconds


%setup basic parameters for the flow
cso  = [1,1,1];
csi  = [64,64,1];
F    = 50/100/60*hd; %absolute in-/outflow in mm^3/s
K    = 5e-6;         %Units in mm^2
mu   = 5*1e-3;       %viscosity is approx 5 cP = 5e-3 Poise


%I tried to calculate the permeability with the formula: K=(1/24)*pi*n*r^4;
%Here n is the number of vessels?
%Following a source (see publication) the av. capillary radius is about 3microns
%{
CBV  = 0.05; %in percent
r0   = 3e-3; %in mm
vCap = r0^2*pi*h(1); %average volume of capillary vessel
n    = hd*CBV/vCap;
K    = 1/24*pi*n*r0^4;
%}


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


%% do the job
if renewCalculation
    
    
    %% setup the flow

    %in-/outflow
    FMat          = zeros(m); 
    Fmat(1,1,1)   = F; 
    Fmat(64,64,1) = -F;

    %permeability
    Kmat          = K*ones([m,3]); %Kmat is a (diagonal) tensor

    %get the ultimate flow
    [pmat,qmat]      = perfusion1c.syntTPFA(Kmat,Fmat,prm); %qmat is absolute flow, i.e. qmat = -K/mu * \nabla p * dA

    %% setup the transport

    %setup AIF
    switch aif
        case 'delta'
            aifval    = zeros(size(timeline));
            aifval(1) = phimat(1);

        case 'gamma'
            aifval = perfusion1c.getGammaAIF(timeline/60); %the AIF is assumed to be in mMol/l.
            aifval = aifval*1e-6;                          %Conversion to mMol/mm^3 by multiplication with 1e-6

    end
    


    %run forward problem    
    disp('Starting simulation...')
    Cmat = perfusion1c.syntforwprobpde(phimat,qmat,Fmat,aifval,prm);
    disp('....done!');



    %downsample the data to approx 1sec time resolution
    step     = round(1/dt);
    Cmat     = Cmat(:,:,:,1:step:end);
    timeline = timeline(1:step:end);
    k        = numel(timeline);


end


%% show the flow



if dynamicFlow

    %display the flow
    cmax = .1*max( Cmat(:) );

    figure(1);clf;
    colormap gray;
    set(1,'name','Idicator-Flow comparison')
    
    for i = 1:k

        Ci  = squeeze(Cmat(:,:,1,i));

        imagesc(Ci);
        axis image
        caxis([0,cmax]);
        ti = sprintf('PDE-Indicator, %1.0f%% completed',i/k*100);
        title(ti);
        
        drawnow;

    end

end

