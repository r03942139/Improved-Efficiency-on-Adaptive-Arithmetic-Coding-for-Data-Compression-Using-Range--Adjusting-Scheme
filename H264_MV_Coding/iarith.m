function x = iarith(code, P, N, sign, method)
    top0 = 1;
    bot0 = 0;
%     P=ones(1,R(2));
    C = [0 cumsum(P)]/sum(P);
    cod = code-48;         % modification
    tag1 = 0.5*cod(1);     % modification
    tag2 = tag1+0.5;       % modification
    b = 1;           % modification
    k = 1;           % modification
    x = zeros(1,N);
    for i = 1:N
        check = true;
        range = bot0 + (top0-bot0) *C;
        j = 1;
        j1 = 0;
        
        while check
            if (range(j+1)>=tag1) && (j1==0);  % modification
                j1 = j;
            end
            if range(j) <= tag1 && range(j+1) >= tag2  % modification
                x(i) = j;
                bot0 = range(j);
                top0 = range(j+1);
                [C,P] = calProb(P, x(i), sign, method);  % update prob.
                check = false;
            else
                j = j+1;
            end
            if range(j) >= tag2  % modification
                j = j1; % j=j-1;
                k = k+1;
                b = b+1;
                tag1 = tag1+2^(-b)*cod(k);
                tag2 = tag1+2^(-b);
            end
        end   
        
        check0 = top0<=0.5;
        check1 = bot0>=0.5;
        while check0 || check1
            b = b-1;  % modification
            if check0
                top0 = top0*2;
                bot0 = bot0*2;
                tag1 = tag1*2;   % modification
                tag2 = tag2*2;   % modification              
                % out=[out '0'];
            elseif check1
                top0 = (top0-0.5)*2;
                bot0 = (bot0-0.5)*2;
                tag1 = (tag1-0.5)*2; % modification
                tag2 = (tag2-0.5)*2; % modification
                % out=[out '1'];
            end
            check0 = top0<=0.5;
            check1 = bot0>=0.5;
        end
    end
end