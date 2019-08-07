function E=calEntropy(x,R)
    P=zeros(1,R);
    for i=1:R
        P(i)=sum(x==i);
    end
    P=P/length(x);
    
    E=-sum(P.*log2(P+eps));

end