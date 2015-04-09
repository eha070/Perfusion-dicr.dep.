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


dynamicFlow = 0;
showMaps    = 0;
showCurve   = 1;


%% load data
prm = settings;
basename   = perfusion1c.struct2nameIndicator(prm,'phiopt','Kopt','dim','aiftype','T','stepred');
foldername = './results/';

%load PDE indicatorflow, CmatPDE
fname = [foldername,'synt-createindicatorpde-' basename '.mat'];
msg      = ['Loading ' fname];
disp(msg);
load(fname);
CmatPDE = Cmat; clearvars Cmat;
timeline = prm.timeline;


%load convolution indicatorflow CmatConv
fname = [foldername,'syntconv-createindicatorconv-' basename '.mat'];
msg   = ['Loading ' fname];
disp(msg);
load(fname,'Cmat');
CmatConv = Cmat; clearvars Cmat;



%load CBF,CBV
prm = settings;
basenameFlow = perfusion1c.struct2nameIndicator(prm,'phiopt','Kopt','dim');
fname = [foldername,'synt-createflowTPFA-' basenameFlow '.mat'];
msg   = ['Loading ' fname];
disp(msg);
load(fname,'phimat','perfmat','qmat');


%setup CBV
CBV  = phimat;


%setup CBFPerf and CBFav, the "averaged" CBF
CBFPerf  = perfmat;
qcc = perfusion1c.convertFlowStagToCC(qmat);
CBFAv = .5*(qcc{1} + qcc{2});
CBFAv(1,1) = 2*CBFAv(1,1);
CBFAv(end,end) = 2*CBFAv(end,end);

clearvars phimat perfmat qmat;



%get main variables
k    = numel(aifval);

%% dynamic flow

if dynamicFlow

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


if showCurve

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
    
    

    MTTPerf = CBV./CBFPerf;
    MTTAv   = CBV./CBFAv;
    
    %setup scaling
    cmaxMTTPerf = .05*max(MTTPerf(:));
    cmaxMTTAv   = .05*max(MTTAv(:));
    cmaxCBFPerf = .05*max(100*60*CBFPerf(:));
    cmaxCBFAv   = .05*max(100*60*CBFAv(:));
    
    cmaxCBF     = max(cmaxCBFAv,cmaxCBFPerf);
    cmaxMTT     = max(cmaxMTTPerf,cmaxMTTAv);


    figure(3);clf; 
    colormap jet(512);
    set(3,'name','Comparison in Flow and MTT')

    subplot(2,2,1)
    imagesc(MTTPerf);
    axis image;
    caxis([0,cmaxMTT])
    title('MTT (perfusion), sec');
    
    subplot(2,2,3);
    imagesc(CBFPerf*100*60)
    title('Flow (perfusion), mmol/min/100g');
    caxis([0,cmaxCBF])
    axis image;        

    
    subplot(2,2,2)
    imagesc(MTTAv);
    axis image;
    caxis([0,cmaxMTT])
    title('MTT (averaging), sec');    

    subplot(2,2,4);
    imagesc(CBFAv*100*60)
    title('Flow (averaging), mmol/min/100g');
    caxis([0,cmaxCBF])
    axis image;
    

    
    
end



