function [C,P] = calProb(P, a, sign, method)
% ===== Original SAC =====
if strcmp(method, 'SAC') == 1
    C = [0 cumsum(P)]/sum(P);   
    return;
% ===== Original AAC =====
elseif strcmp(method, 'AAC') == 1
    P(a) = P(a)+1;
    C = [0 cumsum(P)]/sum(P); 
% ===== Proposed AAC =====    
elseif strcmp(method, 'Pro') == 1
    P(a) = P(a)+1;
%% ***** Section II. E. Local Frequency Table *****
    if  a <= 5 
        zo = ones(1,length(P));
        if sign == 'x'
            if a == 1
                zo(1) = 6; 
            elseif a == 2
                zo(1) = 2;
                zo(2) = 4;
            elseif a == 3
                zo(1) = 2;
                zo(3) = 3;
            elseif a == 4
                zo(4) = 4;% 4
            elseif a == 5
                zo(5) = 4;
            end
            
        elseif sign == 'y'
            if a == 1
            	zo(1) = 4; 
            elseif a == 2
            	zo(2) = 8;
            	zo(4) = 2;
            elseif a == 3
            	zo(3) = 4;
            elseif a == 4
            	zo(2) = 2;
            	zo(4) = 4;
            elseif a == 5
                zo(5) = 4;
            end
        end
        
        PL = P.*zo;        
        C = [0 cumsum(PL)]/sum(PL);   
    else
        C = [0 cumsum(P)]/sum(P);        
    end                
end

end