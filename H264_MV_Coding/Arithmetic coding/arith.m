function out=arith(x,P,R)
    top=1;
    bot=0;

    C=[0 cumsum(P)]/sum(P);
    out=[];
    
    for i=1:length(x)
        [top,bot]=divRange(top,bot,C,x(i));
        check0=top<=0.5;
        check1=bot>=0.5;
        while check0 || check1
            if check0
                top=top*2;
                bot=bot*2;
                out=[out '0'];
            elseif check1
                top=(top-0.5)*2;
                bot=(bot-0.5)*2;
                out=[out '1'];
            end
            check0=top<=0.5;
            check1=bot>=0.5;
        end
        [C,P]=calProb(P,x(i),R);  % update prob.
    end
    ck=1;                 %% modification
    b1=0;
    while ck
        b1=b1+1;
        b2=2^b1;
        c1=ceil(bot*b2);
        c2=floor(top*b2);
        ck=(c1==c2);
    end
    out=[out dec2bin(c1,b1)];
end
