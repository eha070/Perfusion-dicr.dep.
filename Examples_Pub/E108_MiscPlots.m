% ------------------------------------------------------------------------- 
%                               E108_MiscPlots 
% 
% Miscellaneous plots
%
% 1) Quiver-Plot for qmat
% 2) Image of perfmat
% 3) Surface of pressure
% 
%                                          (c)Constantin Heck, 23-Feb-2015 
%                                                http://mic.uni-luebeck.de
% ------------------------------------------------------------------------- 




clc;
close all;
clear;

plotqmat = 0;
plotPerfmat = 0;
plotLperfmat = 1;


%% setup data

load smallDataSet.mat
h = prm.h;


%convert ot cc
qcc = perfusion1c.convertFlowStagToCC(qmat);
qy  = -qcc{1};
qx  = qcc{2};

%setup domain
hx   = h(2);
hy   = h(1);

miny = hy/2;
minx = hx/2;
maxy = m(1)*hy-hy/2;
maxx = m(2)*hx-hx/2;

ygv   = fliplr(miny:hy:maxy);
xgv   = (minx:hy:maxx);
[X,Y] = meshgrid(xgv,ygv);



%% quiver plot of qmat

if plotqmat


    %Quiver plot
    fig = figure(1);clf;
    quiver(X(:),Y(:),qx(:),qy(:),10)
    xlabel('mm');
    ylabel('mm');
    set(gca,'FontSize',15);
    axis image
    xlim([0,3]);
    ylim([0,3]);
    
    export_fig qmat.pdf -transparent

end



%% display perfusion

if plotPerfmat
    
    figure(2);clf;
    imagesc(xgv,ygv,perfmat*100*60);
    axis image;
    xlabel('mm');
    ylabel('mm');    
    set(gca,'YDir','normal','FontSize',15)
    colorbar('south','color','black')

    export_fig perfmat.pdf -transparent

end




%% display local

if plotLperfmat
    
    ca = 100*60*[min(lperfmat(:)),.1*max(lperfmat(:))];
    
    figure(4);clf;
    imagesc(xgv,ygv,lperfmat*100*60);
    caxis(ca);
    axis image;
    xlabel('mm');
    ylabel('mm');    
    set(gca,'YDir','normal','FontSize',15)
    colorbar('south','color','black')

    export_fig lperfmat.pdf -transparent

end
