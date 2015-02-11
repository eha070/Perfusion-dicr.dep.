function [] = cell2tex(A,filename,n)
% CELL2TEX prints a cell array to a latex table
% 
% CELL2TEX(A,FILENAME,N) prints the cell array A to a latex table in
% textfile filename with n digits. The input is a cell array to be able to print
% strings also
%
% Ex: table like
% \begin{center}
%   \begin{tabular}{ l | c || r | }
%     \hline
%     1 & 2 & 3 \\ \hline
%     4 & 5 & 6 \\ \hline
%     7 & 8 & 9 \\
%     \hline
%   \end{tabular}
% \end{center}
%
%

fid = fopen(filename,'wt');

[M N] = size(A);
str = '\begin{table}[H]';
fprintf(fid,'%s\n',str);
str = '\scriptsize';
fprintf(fid,'%s\n',str);
str = '\begin{center}';
fprintf(fid,'%s\n',str);
e = [];
for i = 1 : size(A,2)
    e = [e ' l '];
end;
str = ['\begin{tabular}{' e '}'];
fprintf(fid,'%s\n',str);
format = ['%6.' int2str(n) 'f'];
for i = 1 : M
    for j = 1 : N
        if isnumeric(A{i,j})
            fprintf(fid,format,A{i,j});
        else
            a = A{i,j};
            ind = findstr(a,'_');
            a(ind) = ' ';
            fprintf(fid,'%s',a);
        end;
        if j < N
            fprintf(fid,'%s',' & ');
        end;
    end;
    fprintf(fid,'%s\n','\\');
    if i == 1
        str = '\hline';
        fprintf(fid,'%s\n',str);
    end;
end;
str = '\hline';
fprintf(fid,'%s\n',str);

[M N] = size(A);
str = '\end{tabular}';
fprintf(fid,'%s\n',str);
str = '\end{center}';
fprintf(fid,'%s\n',str);
str = '\caption{}';
fprintf(fid,'%s\n',str);
str = '\label{}';
fprintf(fid,'%s\n',str);
str = '\end{table}';
fprintf(fid,'%s\n',str);

fclose(fid);

