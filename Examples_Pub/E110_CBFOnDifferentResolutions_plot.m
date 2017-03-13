function [] = E110_CBFOnDifferentResolutions_plot

pathload = 'E110_CBFOnDifferentResolutions-results.mat';
D = load(pathload);
results = D.results;
h = results.prm.h;

% Convert to ml/min/100ml
scaling = 60*100;

n = numel(results.P);
for i = 1 : n
    P = results.P{i};
    PLocal = results.PLocal{i};
    CBFCirc = results.CBFCirc{i};
    CBFMS = results.CBFMS{i};
    
    % Raw perfusion 
    val = D.results.P{i}*scaling;
    perf.P(i,1) = mean(val(:));
    
    val = D.results.PLocal{i}*scaling;
    perf.PLocal(i,1) = mean(val(:));
    
    val = D.results.CBFCirc{i}*scaling;
    perf.circ(i,1) = mean(val(:));
    
    val = D.results.CBFMS{i}*scaling;
    perf.ms(i,1) = mean(val(:));
    
    % For CBFCirc
    a = 100*(CBFCirc - P)./CBFCirc;
    err.circ.P.relerror(i,1) = mean(a(:));          
    err.circ.P.absrelerror(i,1) = mean(abs(a(:)));
    
    a = 100*(CBFCirc - PLocal)./CBFCirc;
    err.circ.PLocal.relerror(i,1) = mean(a(:));          
    err.circ.PLocal.absrelerror(i,1) = mean(abs(a(:)));

    % For CBFMS
    a = 100*(CBFMS - P)./CBFMS;
    err.ms.P.relerror(i,1) = mean(a(:));          
    err.ms.P.absrelerror(i,1) = mean(abs(a(:)));

    a = 100*(CBFMS - PLocal)./CBFMS;
    err.ms.PLocal.relerror(i,1) = mean(a(:));          
    err.ms.PLocal.absrelerror(i,1) = mean(abs(a(:)));

    
end;

% Plot errors
H = h(1)*results.voxelSizeFactorList;

symb = {':sk',':ob',':db','-or','-dr'};
figsize = [500,400,500,300];

%
% Plot for Ps
%
handlefig = figure('Position',figsize);

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

xlabel('Voxel size (mm)');
ylabel('Perfusion (ml/min/100ml)')

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

xlabel('Voxel size (mm)');
ylabel('Average reconstruction error RE (%)')

legend('P_s','P_{bSVD}','P_{MS}','RE(bSVD)','RE(MS)')
hold off;

pathsave = ['../pub/figs/' mfilename '-Ps.eps'];
msg = ['Saving ' pathsave];
disp(msg);
print(pathsave,'-depsc');

%
% Plot for Pv
%
handlefig = figure('Position',figsize);

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
xlabel('Voxel size (mm)');
ylabel('Logarithm of perfusion (ml/min/100ml)')

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
xlabel('Voxel size (mm)');
ylabel('Logarithm of average reconstruction error RE (%)')

legend('P_v','P_{bSVD}','P_{MS}','RE(bSVD)','RE(MS)')
hold off;

pathsave = ['../pub/figs/' mfilename '-Pv.eps'];
msg = ['Saving ' pathsave];
disp(msg);
print(pathsave,'-depsc');
