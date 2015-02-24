% Made by Constantin Heck
%
% Modified by Erlend Hodneland
%

function A = getCircularConvolutionMatrix(AIF,deltaT)


    k = numel(AIF);
    AIFHat = [AIF(:);zeros(k,1)];
    
    %setup main variables
    n = numel(AIFHat);

    %standard blockcirculant matrix
    A  = AIFHat;
    for i = 1:(n-1)
        A      = [A,circshift(AIFHat,i)];
    end

%     %Implementation from Wu et al.
%     A      = zeros(n,n);
%     for i = 1:n
%         for j = 1:n
%             if j<=i
%                 A(i,j) = AIFHat(i-j+1);
%             else
%                 A(i,j) = AIFHat(n+(i-j));
%             end
%         end
%     end
    
    %include integration
    A = deltaT*A;



end
