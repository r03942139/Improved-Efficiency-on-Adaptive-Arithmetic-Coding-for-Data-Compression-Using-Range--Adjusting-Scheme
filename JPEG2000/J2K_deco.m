function out = CXdeco(code,CX)
%%  Context table setup
history=12;
k=history+1;
counts=ones(19*k,2);
for i=1:k
    counts(i:19:end,1)=counts(i:19:end,1)+k-i;
end
g = fspecial('gaussian',[1,2*k],2.9);
g = g(k+1:end)/max(g);
%% Parameters
data_high = 1; %the temp of upper bound
data_low  = 0; %the temp of lower bound
bit_high  = 1; %the second temp of upper bound
bit_low   = 0; %the second temp of lower bound
N=length(CX);
out=zeros(2,N);
code_loc=0; % record the location of code
bit_scale=0.5; % scaling number
ADD=1;
%% Table switching parameters
CX_index=CX+1;
table_index=1;
%% Main
for i = 1:N
    context=CX_index(i);
    flag=1;
    while(flag)
        s=cumsum(counts(context*k-k+table_index,:)/sum(counts(context*k-k+table_index,:)));
        prevupper=data_high;
        prevlower=data_low;
        for j=1:size(counts,2)% j=1,2
            templower=data_low;
            if j~=1 % j=2
                templower=templower+(prevupper-prevlower)*s(j-1);
            end
            tempupper=prevlower+(prevupper-prevlower)*s(j);
            if templower<=bit_low && bit_high<=tempupper
                out(2,i)=j;
                flag=0;
                data_high=tempupper;
                data_low=templower;
            end
        end
        if(flag == 0)
            break;
        end
        code_loc=code_loc+1;
        if code_loc>length(code)
            code(code_loc)=1;
        end
        if code(code_loc)==1
            bit_low  = bit_low+bit_scale;
        else
            bit_high = bit_high-bit_scale;
        end
        bit_scale=bit_scale*0.5;
    end
%% Range adjusting
    for j=1:k
        dist=abs(table_index-j)+1;
        counts(context*k-k+j,out(2,i))=counts(context*k-k+j,out(2,i))+1.7*ADD*g(dist);
    end
%% Mutual learning
    % Only for Zero Coding Contexts, and only learns from neighboring contexts
    if context~=1 && context<=9 
        for j=1:k
            dist=abs(table_index-j)+1;
            counts((context-1)*k-k+j,out(2,i))=counts((context-1)*k-k+j,out(2,i))+ADD*g(dist)/5;
        end
    end
    if context~=19 && context<=9
        for j=1:k
            dist=abs(table_index-j)+1;
            counts((context+1)*k-k+j,out(2,i))=counts((context+1)*k-k+j,out(2,i))+ADD*g(dist)/5;
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
                data_high=2*data_high;
                data_low=2*data_low;
                bit_high=2*bit_high;
                bit_low=2*bit_low;
          else
                data_high = 2*data_high -1;
                data_low  = 2*data_low -1;
                bit_high  = 2*bit_high-1;
                bit_low   = 2*bit_low-1;
          end
          bit_scale = 2*bit_scale;   %bit backward
          if data_high == data_low
              fprintf('ERROR: upper = lower occurs !\n');
              return;
          end
    end
%% Table switching
    if i < history
        num = sum(out(2,1:i)==1);
    else
        num = sum(out(2,i-history+1:i)==1);
    end
    table_index = num+1;
end
out=out-1;
out(1,:)=CX(1,:);

end

