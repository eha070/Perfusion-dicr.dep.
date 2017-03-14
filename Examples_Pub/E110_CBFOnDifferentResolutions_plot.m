function [] = E110_CBFOnDifferentResolutions_plot(scaleto,islegend)

pathload = 'E110_CBFOnDifferentResolutions-results.mat';
D = load(pathload);
results = D.results;
h = results.prm.h;
voxelSizeFactorList = D.results.voxelSizeFactorList;

% % Convert to ml/min/100ml
% scaling = 60*100;

% Scale to volume (V), surface (S) or none?
% scaleto = 'none';
% scaleto = 'V';
% scaleto = 'S';

n = numel(results.P);
for i = 1 : n
    P = results.P{i};
    PLocal = results.PLocal{i};
    CBFCirc = results.CBFCirc{i};
    CBFMS = results.CBFMS{i};
    
    % Voxel size
    H = voxelSizeFactorList(i)*h;
    H(3) = h(3);
        
    % Volume in mm^3
    V = prod(H);
    
    % Surface in mm^2
    S = H(1)*H(2)*2 + H(1)*H(3)*2 + H(2)*H(3)*2;

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
        
        ylab1 = 'Q/S (ml/min/mm^2)';                
        
    elseif isequal(scaleto,'none')
        
        %
        % Convert from mm^3/s/mm^3 to ml/min/100ml        
        %
        
        % mm^3/s/mm^3 = ml/s/ml, nothing to do        
        
        % From ml/s/ml to ml/min/ml
        scaling = 60;
        
        % From ml/min/ml to ml/min/100ml
        scaling = scaling*100;
        
        ylab1 = 'Q/100ml (ml/min/100ml)';
    else
        error('Use option none or S for scaling');
    end;
    
    %
    % Raw perfusion 
    %
    val = D.results.P{i}*scaling;    
    perf.P(i,1) = mean(val(:));    
    
    val = D.results.PLocal{i}*scaling;
    perf.PLocal(i,1) = mean(val(:));
    
    val = D.results.CBFCirc{i}*scaling;    
    perf.circ(i,1) = mean(val(:));
    
    val = D.results.CBFMS{i}*scaling;    
    perf.ms(i,1) = mean(val(:));
    
    %
    % Error
    %
    
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

xlabel('Voxel size (mm)');
ylabel(ylab1)

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
xlabel('Voxel size (mm)');
ylabel(['Logarithm of ' ylab1]);

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

if islegend
    legend('P_v','P_{bSVD}','P_{MS}','RE(bSVD)','RE(MS)')
end;
hold off;

pathsave = ['../pub/figs/' mfilename '-Pv' '-scaleto-' scaleto '.eps'];
msg = ['Saving ' pathsave];
disp(msg);
print(pathsave,'-depsc');
