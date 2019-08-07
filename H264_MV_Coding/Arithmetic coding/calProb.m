function [C,P]=calProb(P,a,flag)
% flag: # of data


% if a<10
    P(a) = P(a)+1;
% else
%     P(a) = P(a)+4;
% end
    
 
    C=[0 cumsum(P)]/sum(P);
end