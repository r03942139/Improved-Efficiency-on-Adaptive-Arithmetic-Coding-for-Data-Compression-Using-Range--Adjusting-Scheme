function out = arith_Wx(x,P)
    top = 1;
    bot = 0;

    C = [0 cumsum(P)]/sum(P);
    out=[];
    
    for i = 1:length(x)
%%        [top,bot]=divRange(top,bot,C,x(i));
        bot = bot+(top-bot)*C(x(i));
        top = bot+(top-bot)*C(x(i)+1);
        check0 = top<=0.5;
        check1 = bot>=0.5;
        while check0 || check1
            if check0
                top = top*2;
                bot = bot*2;
                out = [out '0'];
            elseif check1
                top = (top-0.5)*2;
                bot = (bot-0.5)*2;
                out = [out '1'];
            end
            check0 = top<=0.5;
            check1 = bot>=0.5;
        end
        
%%         [C,P] = calProb(P,x(i),R);  % update prob.
%% Increasing Step Size
            if x(i) <= 1 && i >= length(x)*0.69 
                A = 1 + round(0.007*sqrt(i));
            else
                A = 1;
            end
            P(x(i)) = P(x(i)) + A;

%% Rear symbol tunung 
%{
            len = length(P);
            for kk = 1:len
                as = sum(P(1:kk));
                bs = as + sum(P(kk+1:end));
%                 Pro = sum(P(1:kk))/sum(P);
                Pro = as/bs;
                if Pro >= 0.79
                    MSS = kk;
                    break;
                end
            end
            
            if i == round(length(x)/2) 
%                 P(1:10) = round(P(1:10)*3/4)+1;
                P(1:MSS) = round(P(1:MSS)*3/4)+1;
%             elseif i == round(length(x)*3/4)
%                 P = round(P/2)+1;
            end
%}            
%% Local Frequency Table

            if  x(i) <= 5
                zo = ones(1,length(P));
                if x(i) == 1
                    zo(1) = 6; 
                elseif x(i) == 2
                    zo(1) = 2;
                    zo(2) = 4;
                elseif x(i) == 3
                    zo(1) = 2;
                    zo(3) = 3;
                elseif x(i) == 4
                    zo(4) = 4;% 4
                elseif x(i) == 5
                    zo(5) = 4;
                end
                PL = P.*zo;        
                C = [0 cumsum(PL)]/sum(PL);               
            else

                C = [0 cumsum(P)]/sum(P);        
            end            
    end
    ck = 1;                 %% modification
    b1 = 0;
    while ck
        b1 = b1+1;
        b2 = 2^b1;
        c1 = ceil(bot*b2);
        c2 = floor(top*b2);
        ck = (c1==c2);
    end
    out = [out dec2bin(c1,b1)];
    %figure; plot(P)    
end


    
