function [out, x_end] = arith(x, P, sign, method)
    top=1;
    bot=0;
%     P=ones(1,R(2));
    C=[0 cumsum(P)]/sum(P);
    out=[];
    
    for i = 1:length(x)
        [top,bot] = divRange(top,bot,C,x(i));
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
        if strcmp(method, 'SAC') == 0
            [C,P] = calProb(P, x(i), sign, method);  % update prob.
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
    x_end = length(out);
end
