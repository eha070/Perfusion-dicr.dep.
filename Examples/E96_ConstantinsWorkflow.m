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
dynamicFlow      = 1;
setupPDE         = 1;
setupConvolution = 0;
downsampling     = 1;
makelenmat       = 1;
streamlines      = 'lenmatE';
saveDataSmall    = 1;

%configuration
% aif = 'delta';
aif = 'gamma';
% aif = 'parker';
% aif = 'test';
% aif = 'delta2';


%setup spatial-temporal domain
m        = [64,64,1];       %matrix size
omega    = [0,10,0,10,0,1]; %domain in mm
Tmax     = 90;              %max time in seconds
dt       = 0.0005;          %time resolution of simulation
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
K    = 5e-6;         %Units in mm^2
mu   = 5*1e-3;       %viscosity is approx 5 cP = 5e-3 Pa*s


%I tried to calculate the permeability with the formula: K=(1/24)*pi*n*r^4;
%Here n is the number of vessels?
%Following a source (see publication) the av. capillary radius is about 3microns
%{
r0   = 3e-3;         %in mm
vCap = r0^2*pi*h(1); %average volume of capillary vessel
CBV  = 0.05;         %in percent
n0   = hd*CBV/vCap;
K    = 1/24*pi*n0*r0^4;
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
    [pmat,qmat]      = perfusion1c.syntTPFA(Kmat,Fmat,prm); %qmat is absolute flow, i.e. qmat = -K/mu * \nabla p * dA

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


%% create the streamlines

if makelenmat

    switch streamlines

        case 'lenmatC'
            %flow: Go from stag to cc
            qcc = perfusion1c.convertFlowStagToCC(qmat);
            qy  = qcc{1};
            qx  = qcc{2};

            %get meshgrid to plot qcc, the cell-centered flow
            ygv   = linspace(omega(2),omega(1),m(1)); %from (0 to 1)
            xgv   = linspace(omega(3),omega(4),m(2)); %this goes from (1 to 0) (matlab, duh)
            [x,y] = meshgrid(xgv,ygv);

            %seed-points
            coord = zeros(n,2);
            for i = 1:m(1);
                for j = 1:m(2);
                    idx          = sub2ind(m,i,j);
                    coord(idx,:) = [h(1)/2 + (i-1)*h(1), h(2)/2 + (j-1)*h(2)];
                end
            end

            %get the streamlines
            quiver(x(:),y(:),qx(:),-qy(:),10)
            fh = streamline(x,y,qx,-qy,coord(:,2),coord(:,1),[.01,100000]); %streamlines upstream

            %get the lengths of the streamlines
            l = zeros(n,1);
            for i = 1:n;

                %get upstream streamline length
                resx = fh(i).XData(2:end)-fh(i).XData(1:end-1);
                resy = fh(i).YData(2:end)-fh(i).YData(1:end-1);
                l(i) = sum(sqrt(resx.^2+resy.^2));

            end
            lenmat = reshape(l,m);
            lenmat = flipud(lenmat);
        case 'lenmatE'
            lenmat = perfusion1c.arclength(qmat,Fmat,h);
    end
    
    
    
else
    
    load('tmp.mat','lenmat');

    
end
    



% get perfusion

%get norm of NORMALIZED flow. 
qcc  = perfusion1c.convertFlowStagToCC(qmat);
A1 = h(2)*h(3);
A2 = h(1)*h(3);    

qy   = qcc{1}./A1;
qx   = qcc{2}./A2;
qabs = sqrt(qx.^2 + qy.^2);

%get perfusion
CBF = qabs./lenmat;

%get CBV
CBV = phimat;
            
 

%% save some data
if saveDataSmall
    
    %downsample the data to approx 1sec time resolution
    step     = round(tSamp/dt);
    Cmat     = Cmat(:,:,:,1:step:end);
    timeline = timeline(1:step:end);
    aifval   = aifval(1:step:end);
    k        = numel(timeline);
    
    save('smallDataSet.mat','Cmat','timeline','aifval','k','n','Hd','qmat');
    return;
    
end




%% create convolution model

if setupConvolution
    
    
    %setup MTT
    MTT = CBV./CBF;
    

    %setup delay map
    [~,delay]    = max(Cmat,[],4);
    [~,delayaif] = max(aifval);
    delay        = delay - delayaif;

    %shift the AIF
    aifmat = zeros(mk);
    for i = 1 : m(1)
        for j = 1 : m(2)
            for l = 1 : m(3)
                t0 = delay(i,j,l);
                nAIF = k-t0;
                aifmat(i,j,l,t0+1:end) = aifval(1:nAIF);
                aifmat(i,j,l,1:t0) = 0;
            end;
        end;
    end;

    %do the convolution
    dt    = timeline(2)-timeline(1);
    CConv = zeros(size(Cmat));
    
    fprintf('i=');
    for i = 1 : m(1)
        fprintf('%i, ',i);
        for j = 1 : m(2)
            for l = 1 : m(3)
                val = squeeze(aifmat(i,j,l,:));
                IH  = exp(-timeline./MTT(i,j,l));
                Ci   = conv(IH,val);
                Ci   = CBF(i,j,l)*dt*Ci(1:k);       
                CConv(i,j,l,:) = Ci;
            end
        end
    end
    fprintf('\n...done\n');


end


%% downsample the data

if downsampling

    %downsample the data to approx 1sec time resolution
    step     = round(tSamp/dt);
    Cmat     = Cmat(:,:,:,1:step:end);
    CConv    = CConv(:,:,:,1:step:end);
    timeline = timeline(1:step:end);
    k        = numel(timeline);
    aifval   = aifval(1:step:end);
    mk       = [m,k];
    
end



 %% show the flow

if dynamicFlow

    %display the two indicator flows
    cmax = .1*max( [Cmat(:);CConv(:);] );

    figure(1);clf;
    colormap gray;
    set(1,'name','Idicator-Flow comparison')
    
    for i = 1:k

        CiPDE  = squeeze(Cmat(:,:,1,i));
        CiConv = squeeze(CConv(:,:,1,i));

        subplot(1,2,1);
        imagesc(CiPDE);
        axis image
        caxis([0,cmax]);
        ti = sprintf('PDE-Indicator, %1.0f%% completed',i/k*100);
        title(ti);

        subplot(1,2,2);
        imagesc(CiConv);
        axis image
        caxis([0,cmax]);
        ti = sprintf('CONV-Indicator using Perfusion, %1.0f%% completed',i/k*100);
        title(ti);    
        
        
        pause(.1);

        drawnow;

    end


end


