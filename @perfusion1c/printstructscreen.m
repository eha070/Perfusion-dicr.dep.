function [] = printstructscreen(var)
% PRINTSTRUCTSCREEN prints a struct array to the screen
%
%    PRINTSTRUCTSCREEN(VAR) prints the struct array in VAR to the screen
%
%
p = 20;
f = fieldnames(var);
n = 0;
for i = 1 : length(f)
    n = max(n,length(f{i}));    
end;

for i = 1 : length(f)
    val = var.(f{i});
    if numel(val) > 20 && (isnumeric(val) || islogical(val))
        val = [num2str(val(1)) ' to ' num2str(val(end))];
    end;
    if isnumeric(val) || islogical(val);
        val = num2str(val);
    end;

    if iscell(val) || isstruct(val)
        continue;
    end;
    try
        msg = [makestr(f{i},n+3) makestr(':',n+3) makestr(val,n+3)];
        disp(msg);
    catch
        
    end;
    
end;