% ------------------------------------------------------------------------- 
%                         E01_CompareIndicatorFlow
% 
% Displays the indicator flow from the PDE and Convolution for comparison
% 
% 
% 
%                                          (c)Constantin Heck, 30-Jan-2015 
%                                                http://mic.uni-luebeck.de
% ------------------------------------------------------------------------- 

clear;
clc;
close all;


showFlow   = 1;
showMaps   = 0;
showCurves = 0;


%% load data
prm = settings;
basename = perfusion1c.struct2nameIndicator(prm,'phiopt','Kopt','dim','aiftype','T','stepred');

%load PDE indicatorflow
fname = ['synt-createindicatorpde-' basename '.mat'];
msg      = ['Loading ' fname];
disp(msg);
load(fname);
CmatPDE = Cmat; clearvars Cmat;
timeline = prm.timeline;


%load convolution indicatorflow using Erlend's perfusion
fname = ['syntconv-createindicatorconv-' basename '.mat'];

msg      = ['Loading ' fname];
disp(msg);
load(fname,'Cmat');
CmatConv = Cmat; clearvars Cmat;



%load CBV (i.e. phimat)
prm = settings;
basenameFlow = perfusion1c.struct2nameIndicator(prm,'phiopt','Kopt','dim');
fname = ['synt-createflowTPFA-' basenameFlow '.mat'];

msg   = ['Loading ' fname];
disp(msg);
load(fname,'phimat','perfmat');
CBF  = perfmat;
CBV  = phimat; clearvars phimat perfmat;



%get main variables
k    = numel(aifval);



if showFlow

    %display the two indicator flows
    step = 1;
    cmax = .1*max( [CmatConv(:);CmatPDE(:);] );

    figure(1);clf;
    colormap gray;
    set(1,'name','Idicator-Flow comparison')
    
    for i = 1:step:k

        CiPDE  = squeeze(CmatPDE(:,:,1,i));
        CiConv = squeeze(CmatConv(:,:,1,i));

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

        drawnow;

    end

end

%% comparison on single curves


if showCurves

    i = randi(prm.dim(1)/2);
    j = randi(prm.dim(2)/2);


    CPDEij  = squeeze(CmatPDE(i,j,:));
    CConvij = squeeze(CmatConv(i,j,:));

    scalePerf = max(CConvij)./max(CPDEij);


    figure(2);clf;
    set(2,'name','Comparison of a single tissue-curve')
    
    subplot(1,2,1);
    plot(timeline,CPDEij,timeline,CConvij,timeline,CConvij/scalePerf,'linewidth',2);
    legend('CPDE','CConv','Scaled CPerfusion');
    title(sprintf('Scaled: (i,j)=(%i,%i), fac=%1.4f',i,j,scalePerf));


end

%% show maps with CBF,CBV and MTT

if showMaps
    
    

    MTTConv = CBV./CBF;
    
    %setup scaling
    cmaxMTT = .05*max(MTTConv(:));
    cmaxCBF = .05*max(100*60*CBF(:));


    figure(3);clf; 
    colormap jet(512);
    set(3,'name','Comparison in Flow and MTT')

    subplot(1,2,1)
    imagesc(MTTConv);
    axis image;
    caxis([0,cmaxMTT])
    title('MTT (perfusion), sec');


    subplot(1,2,2);
    imagesc(CBF*100*60)
    title('Flow (perfusion), mmol/min/100g');
    caxis([0,cmaxCBF])
    axis image;
    
    
end



