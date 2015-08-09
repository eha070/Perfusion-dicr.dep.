% ------------------------------------------------------------------------- 
%                               ConvertTimelineToTextfile
% 
% 
% 
% 
% 
%                                          (c)Constantin Heck, 06-Aug-2015 
%                                                http://mic.uni-luebeck.de
% ------------------------------------------------------------------------- 

clear;
clc;

load syntconv-createindicatorconv-phi-flat-K-flat-dim-64-aif-gamma-T-90-red-100;



%% write timeline
%initialize variables
timeline    = prm.timeline;
k           = numel(timeline);
relTimeline = timeline(2:end)-timeline(1); 

%initialize variables
s = '2015-01-01 00:00:00.000';


%create string
for i = 1:k-1

    stringIncrement = sprintf(' / +%1.3f',relTimeline(i));
    s               = [s,stringIncrement];
    
end


%save to text file
f = fopen('timeline.txt','w');
fprintf(f,'%s',s);
fclose(f);



%% write aifval
timeline    = prm.timeline;
k           = numel(timeline);
relTimeline = timeline(1:end)-timeline(1); 

s = sprintf('#Measured AIF\n');

%create string
for i = 1:k

    stringIncrement = sprintf('%1.3f    %1.14f\n',relTimeline(i),aifval(i));
    s               = [s,stringIncrement];
    
end


%save to text file
f = fopen('aif.txt','w');
fprintf(f,'%s',s);
fclose(f);