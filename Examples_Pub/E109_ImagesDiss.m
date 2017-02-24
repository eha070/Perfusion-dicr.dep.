% ------------------------------------------------------------------------- 
%                            E109_ImagesDissertation
% 
% 
% 
% 
% 
%                                      (c)Constantin Sandmann, 24-Feb-2017
%                                                http://mic.uni-luebeck.de
%-------------------------------------------------------------------------


clc;
clear;
close all;


load smallDataSet.mat


plotMovie    = 0;
plotPressure = 0;
plotAIF      = 0;

%% plot movie of CA propagation

if plotMovie

    numPlots = 8;
    idx      = round(linspace(1,150,numPlots));

    %prepare plot
    figure(1);
    clf;
    colormap gray;
    ca       = [min(Cmat(:)),max(Cmat(:))];

    %prepare domain
    xgv = linspace(0,prm.physdim(2),m(2));
    ygv = linspace(0,prm.physdim(1),m(1));


    %plot
    for i = 1:numPlots
        imagesc(Cmat(:,:,1,idx(i)));
        caxis(ca);
        axis off;
        axis image;
        fname = sprintf('SimTrans_%1.0fs.pdf',timeline(idx(i)));
        export_fig(fname,'-transparent');
    end

end



%% display pressure

if plotPressure
    
    pmat = pmat-min(pmat(:));
    figure(2);clf;
    surf(xgv,ygv,pmat);
    axis equal
    xlabel('mm');
    ylabel('mm');    
    zlabel('Pa');    
    set(gca,'YDir','normal','FontSize',15)
    view(110,24)    
    daspect([1 1 2e2]);

    export_fig pmat.pdf -transparent

end


%% display aif

if plotAIF
    
    figure(3);clf;
    plot(timeline,aifval,'linewidth',3);
    legend('AIF');
    xlabel('time [s]');
    ylabel('concentration [mmol/mm^3]');    
    set(gca,'FontSize',15);    
    export_fig synth_aif.pdf -transparent
    
end