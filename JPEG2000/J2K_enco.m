% Adaptive Arithmetic  Encoding  Function
% notes: 0~8 zero coding   9~13: sign coding  14~16: magnitude refinement coding 17~18: run-length coding
function out = CXenco(in)
%%  Context table setup & Initialization
history=12;
k=history+1;
counts=ones(19*k,2);
for i=1:k
    counts(i:19:end,1)=counts(i:19:end,1)+k-i;
end
g = fspecial('gaussian',[1,2*k],2.9);
g = g(k+1:end)/max(g);
%% Parameters
N=length(in);
out=[];
data_high=1;
data_low=0;      
ADD = 1;
%% Table switching parameters
table_index=1;      
CX_index = in(1,:)+1;
D_index  = in(2,:)+1;
%% Main
for i = 1:N
    context=CX_index(i);   
    s=cumsum(counts(context*k-k+table_index,:)/sum(counts(context*k-k+table_index,:)));
    prevupper=data_high;
    prevlower=data_low;
    if D_index(i)~=1
        data_low=data_low+(prevupper-prevlower)*s(D_index(i)-1);
    end
    data_high=prevlower+(prevupper-prevlower)*s(D_index(i));
    
%% Range adjusting
    for j=1:k
        dist=abs(table_index-j)+1;
        counts(context*k-k+j,D_index(i))=counts(context*k-k+j,D_index(i))+1.7*ADD*g(dist);
    end
%% Mutual learning
    % Only for Zero Coding Contexts, and only learns from neighboring contexts
    if context~=1 && context<=9 
        for j=1:k
            dist=abs(table_index-j)+1;
            counts((context-1)*k-k+j,D_index(i))=counts((context-1)*k-k+j,D_index(i))+ADD*g(dist)/5;
        end
    end
    if context~=19 && context<=9
        for j=1:k
            dist=abs(table_index-j)+1;
            counts((context+1)*k-k+j,D_index(i))=counts((context+1)*k-k+j,D_index(i))+ADD*g(dist)/5;
        end
    end
%% Rescaling 
    ADD = ADD*1.031;
    if ADD>2
        counts(context*k-k+1:context*k,:)=ceil(counts(context*k-k+1:context*k,:)/2);
        ADD=ADD/2;
    end
    while (data_high<=0.5)||(data_low>=0.5)
          if data_high<=0.5
               out = [out 0];
               data_high = 2 * data_high;
               data_low =  2 * data_low;
          else
               out = [out 1];
               data_high = 2 * data_high -1;
               data_low =  2 * data_low -1;
          end
          if data_high == data_low
              fprintf('ERROR: upper = lower occurs !\n');
              return;
          end
    end
%% Table switching

    if i < history
        num = sum(D_index(1:i)==1);
    else
        num = sum(D_index(i-history+1:i)==1);
    end
    table_index = num+1;
    
end
%% binary encoding
out = [out 0 1]; 
end
