function [] = E110_CBFOnDifferentResolutions_plot(scaleto)

% Run as 
% E110_CBFOnDifferentResolutions_plot('none')% for volume scaling
% E110_CBFOnDifferentResolutions_plot('S')% for surface scaling

pathload = 'E110_CBFOnDifferentResolutions-results.mat';
D = load(pathload);
results = D.results;
h = results.prm.h;
voxelSizeFactorList = D.results.voxelSizeFactorList;
islegend = 1;

m = D.results.m;
% % Convert to ml/min/100ml
% scaling = 60*100;

% Scale to volume (V), surface (S) or none?
% scaleto = 'none';
% scaleto = 'V';
% scaleto = 'S';

n = numel(results.P);
for i = 1 : n
    
    % Voxel size
    H = voxelSizeFactorList(i)*h;
    H(3) = h(3);
        
    % Volume in mm^3
    V = prod(H);
    
    % Surface in mm^2
    S = H(1)*H(2)*2 + H(1)*H(3)*2 + H(2)*H(3)*2;

    % setup domain for plotting
    miny = h(2)/2;
    minx = h(1)/2;
    maxy = m(1)*h(2)-h(2)/2;
    maxx = m(2)*h(1)-h(1)/2;
    xgv   = (minx:h(1):maxx);
    ygv   = fliplr(miny:h(2):maxy);


    if isequal(scaleto,'S')
        %
        % Convert from mm^3/s/mm^3 to mL/min/mm^2       
        %
        
        % From mm^3/s/mm^3 to mm^3/s
        scaling = V;
        
        % From mm^3/s to mm^3/min
        scaling = scaling*60;
        
        % From mm^3/min to ml/min
        scaling = scaling*1e-3;

        % From ml/min to ml/min/mm^2
        scaling = scaling/S;        
        
        ylab1 = '$$\overline{F}/S$$ $$(ml/min/mm^2)$$';                
        
    elseif isequal(scaleto,'none')
        
        %
        % Convert from mm^3/s/mm^3 to ml/min/100ml        
        %
        
        % mm^3/s/mm^3 = ml/s/ml, nothing to do        
        
        % From ml/s/ml to ml/min/ml
        scaling = 60;
        
        % From ml/min/ml to ml/min/100ml
        scaling = scaling*100;
        
        ylab1 = '$$\overline{F}/100ml$$ $$(ml/min/100ml)$$';
    else
        error('Use option none or S for scaling');
    end;

    P = results.P{i}*scaling;
    PLocal = results.PLocal{i}*scaling;
    CBFCirc = results.CBFCirc{i}*scaling;
    CBFMS = results.CBFMS{i}*scaling;
        
    %
    % Raw perfusion 
    %
    ca1 = [15,70];
    val = P;
    if isequal(scaleto,'none')
        if i == 1
            figure;
            imagesc(xgv,ygv,val);
            axis image;
            xlabel('$$mm$$','interpreter','latex');
            ylabel('$$mm$$','interpreter','latex');    
            set(gca,'YDir','normal','FontSize',22)
            colorbar;
%             colorbar('south','color','black')        
            caxis(ca1);
%             minval = min(val(:));
%             meanval = mean(val(:));
%             maxval = max(val(:));
%             a = ['Min = ' num2str(minval,'%6.0f')];
%             b = ['Mean = ' num2str(meanval,'%6.0f')];
%             c = ['Max = ' num2str(maxval,'%6.0f')];
%             text(1.8,2.9,a,'FontSize',22,'color','white');
%             text(1.8,2.7,b,'FontSize',22,'color','white');
%             text(1.8,2.5,c,'FontSize',22,'color','white');
            pathsave = ['../pub/figs/' mfilename '-Ps' '-scaleto-' scaleto '-raw.eps'];
            msg = ['Saving ' pathsave];
            disp(msg);
            print(pathsave,'-depsc');

        end;
    end;
    perf.P(i,1) = mean(val(:));    
    
    
    ca2 = [1000,20000];
    val = PLocal;
    if isequal(scaleto,'none')
        if i == 1
            figure;
            imagesc(xgv,ygv,val);
            axis image;
            xlabel('$$mm$$','interpreter','latex');
            ylabel('$$mm$$','interpreter','latex');    
            set(gca,'YDir','normal','FontSize',22)
            colorbar;
%             colorbar('south','color','black')        
            caxis(ca2);
%             minval = min(val(:));
%             meanval = mean(val(:));
%             maxval = max(val(:));
%             a = ['Min = ' num2str(minval,'%6.0f')];
%             b = ['Mean = ' num2str(meanval,'%6.0f')];
%             c = ['Max = ' num2str(maxval,'%6.0f')];
%             text(1.8,2.9,a,'FontSize',22,'color','white');
%             text(1.8,2.7,b,'FontSize',22,'color','white');
%             text(1.8,2.5,c,'FontSize',22,'color','white');
            pathsave = ['../pub/figs/' mfilename '-Pv' '-scaleto-' scaleto '-raw.eps'];
            msg = ['Saving ' pathsave];
            disp(msg);
            print(pathsave,'-depsc');            
        end;
    end;
    perf.PLocal(i,1) = mean(val(:));
    
    ca3 = [200,2000];
    val = CBFCirc;
    if isequal(scaleto,'none')
        if i == 1
            figure;
            imagesc(xgv,ygv,val);
            axis image;
            xlabel('$$mm$$','interpreter','latex');
            ylabel('$$mm$$','interpreter','latex');    
            set(gca,'YDir','normal','FontSize',22)
%             colorbar('south','color','black')        
            colorbar;
            caxis(ca3);
%             minval = min(val(:));
%             meanval = mean(val(:));
%             maxval = max(val(:));
%             a = ['Min = ' num2str(minval,'%6.0f')];
%             b = ['Mean = ' num2str(meanval,'%6.0f')];
%             c = ['Max = ' num2str(maxval,'%6.0f')];
%             text(1.8,2.9,a,'FontSize',22,'color','white');
%             text(1.8,2.7,b,'FontSize',22,'color','white');
%             text(1.8,2.5,c,'FontSize',22,'color','white');
            pathsave = ['../pub/figs/' mfilename '-bSVD' '-scaleto-' scaleto '-raw.eps'];
            msg = ['Saving ' pathsave];
            disp(msg);
            print(pathsave,'-depsc');            
        end;
    end;
    perf.circ(i,1) = mean(val(:));
    
    
    ca4 = [40,120];
    val = CBFMS;
    if isequal(scaleto,'none')
        if i == 1
            figure;
            imagesc(xgv,ygv,val);
            axis image;
            xlabel('$$mm$$','interpreter','latex');
            ylabel('$$mm$$','interpreter','latex');    
            set(gca,'YDir','normal','FontSize',22)
            colorbar;
%             colorbar('south','color','red')        
            caxis(ca4);
%             minval = min(val(:));
%             meanval = mean(val(:));
%             maxval = max(val(:));
%             a = ['Min = ' num2str(minval,'%6.0f')];
%             b = ['Mean = ' num2str(meanval,'%6.0f')];
%             c = ['Max = ' num2str(maxval,'%6.0f')];
%             text(1.8,2.9,a,'FontSize',22,'color','white');
%             text(1.8,2.7,b,'FontSize',22,'color','white');
%             text(1.8,2.5,c,'FontSize',22,'color','white');
            pathsave = ['../pub/figs/' mfilename '-MS' '-scaleto-' scaleto '-raw.eps'];
            msg = ['Saving ' pathsave];
            disp(msg);
            print(pathsave,'-depsc');

        end;
    end;
    perf.ms(i,1) = mean(val(:));    
    
    pause
    
    %
    % Error
    %
    
    % For CBFCirc
    a = 100*(CBFCirc - P)./P;
    err.circ.P.relerror(i,1) = mean(a(:));          
    err.circ.P.absrelerror(i,1) = mean(abs(a(:)));
    
    
    a = 100*(CBFCirc - PLocal)./PLocal;
    err.circ.PLocal.relerror(i,1) = mean(a(:));          
    err.circ.PLocal.absrelerror(i,1) = mean(abs(a(:)));

    % For CBFMS
    a = 100*(CBFMS - P)./P;
    err.ms.P.relerror(i,1) = mean(a(:));          
    err.ms.P.absrelerror(i,1) = mean(abs(a(:)));

    a = 100*(CBFMS - PLocal)./PLocal;
    err.ms.PLocal.relerror(i,1) = mean(a(:));          
    err.ms.PLocal.absrelerror(i,1) = mean(abs(a(:)));
    pause
end;

% Plot errors
H = h(1)*results.voxelSizeFactorList;

symb = {':sk',':ob',':db','-or','-dr'};
figsize = [500,400,500,300];

%
% Plot for Ps
%

name = ['Global perfusion, scaling ' scaleto];
handlefig = figure('Position',figsize,'Name',name);

i = 0;
% Left side
i = i + 1;
yyaxis left
handleplot = plot(H,(perf.P),symb{i});
set(handleplot, 'MarkerFaceColor', get(handleplot, 'Color'));
hold on;

yyaxis left
i = i + 1;
handleplot = plot(H,(perf.circ),symb{i});
set(handleplot, 'MarkerFaceColor', get(handleplot, 'Color'));
hold on;

yyaxis left
i = i + 1;
handleplot = plot(H,(perf.ms),symb{i});
set(handleplot, 'MarkerFaceColor', get(handleplot, 'Color'));
hold on;

xlabel('Voxel size ($$mm$$)','Interpreter','latex');
ylabel(ylab1,'Interpreter','latex')

% Right side
yyaxis right
i = i + 1;
handleplot = plot(H,err.circ.P.absrelerror,symb{i});
set(handleplot, 'MarkerFaceColor', get(handleplot, 'Color'));
hold on;

yyaxis right
i = i + 1;
handleplot = plot(H,err.ms.P.absrelerror,symb{i});
set(handleplot, 'MarkerFaceColor', get(handleplot, 'Color'));
hold on;
ylabel('Average reconstruction error RE (\%)','Interpreter','latex')

if islegend
    legend('P_s','P_{bSVD}','P_{MS}','RE(bSVD)','RE(MS)')
end;
hold off;

pathsave = ['../pub/figs/' mfilename '-Ps' '-scaleto-' scaleto '.eps'];
msg = ['Saving ' pathsave];
disp(msg);
print(pathsave,'-depsc');

%
% Plot for Pv
%
name = ['Local perfusion, scaling ' scaleto];
handlefig = figure('Position',figsize,'Name',name);


i = 0;
% Left side
i = i + 1;
yyaxis left
handleplot = plot(H,(perf.PLocal),symb{i});
set(handleplot, 'MarkerFaceColor', get(handleplot, 'Color'));
hold on;

yyaxis left
i = i + 1;
handleplot = plot(H,(perf.circ),symb{i});
set(handleplot, 'MarkerFaceColor', get(handleplot, 'Color'));
hold on;

yyaxis left
i = i + 1;
handleplot = plot(H,perf.ms,symb{i});
set(handleplot, 'MarkerFaceColor', get(handleplot, 'Color'));
hold on;

set(get(handleplot,'Parent'),'YScale','log');
xlabel('Voxel size (mm)','Interpreter','latex');
ylabel(['Logarithm of ' ylab1],'Interpreter','latex');

% Right side
yyaxis right
i = i + 1;
handleplot = plot(H,(err.circ.PLocal.absrelerror),symb{i});
set(handleplot, 'MarkerFaceColor', get(handleplot, 'Color'));
hold on;

yyaxis right
i = i + 1;
handleplot = plot(H,(err.ms.PLocal.absrelerror),symb{i});
set(handleplot, 'MarkerFaceColor', get(handleplot, 'Color'));
hold on;

set(get(handleplot,'Parent'),'YScale','log');
ylabel('Logarithm of average reconstruction error RE (\%)','Interpreter','latex')

if islegend
    legend('P_v','P_{bSVD}','P_{MS}','RE(bSVD)','RE(MS)')
end;
hold off;

pathsave = ['../pub/figs/' mfilename '-Pv' '-scaleto-' scaleto '.eps'];
msg = ['Saving ' pathsave];
disp(msg);
print(pathsave,'-depsc');
