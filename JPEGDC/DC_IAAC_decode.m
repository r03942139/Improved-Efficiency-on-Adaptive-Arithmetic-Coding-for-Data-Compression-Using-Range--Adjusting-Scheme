function dc2 = DC_IAAC_decode(Code,M,N)
% Standard code using the proposed IAAC algorithm to deencode the DC terms 
M1=M/8;
N1=N/8;
MN=M1*N1;
dc2=zeros(1,MN);
Lw=0; Up=1; Lw2=0;  Up2=1;
Le=length(Code);
tk1=min(4000,Le); tk=tk1;
Code1=Code(1:tk);
[vn,Code1,Lw,Up,Lw2,Up2]=dcd(Code1,Lw,Up,Lw2,Up2,ones(1,129)/129);
dc2(1)=vn-1;
tb2=exp(-0.2*[0:128]);
tb2=tb2/sum(tb2)*1000;
tb2=max(tb2,1);
dc1=dc2(1);
add=80;
for a7=2:MN
    dm=max(128-dc1,dc1);
    dm1=min(128-dc1,dc1);
    tb22=tb2(1:dm+1);
    if mod(a7,N)==1
        tb22=tb22.*(1+log2(1:dm+1)*1);
    end
    tb22(dm1+2:dm+1)=tb22(dm1+2:dm+1)*.8;    
    [d1,Code1,Lw,Up,Lw2,Up2]=dcd(Code1,Lw,Up,Lw2,Up2,tb22);
    d2=d1-1;
    if (d2==0)||(dc1-d2<0)
        dc2(a7)=dc1+d2;
    elseif dc1+d2>128
        dc2(a7)=dc1-d2;
    else
        [s1,Code1,Lw,Up,Lw2,Up2]=dcd(Code1,Lw,Up,Lw2,Up2,[.5,.5]);
        if s1 == 2
            dc2(a7)=dc1+d2;
        else
            dc2(a7)=dc1-d2;
        end
    end 
    if d2==0
        ad1=[add,zeros(1,128)];
    else
        d3=abs(([1:128]+1)/(d2+1)-1);
        ad1=[0,exp(-7*d3.^2)];
        ad1=ad1/sum(ad1)*add;        
    end
    tb2=tb2+ad1;
    dc1=dc2(a7);  
    add=add*1.02;            % adjusting step
    if sum(tb2)>100000
        tb2=max(tb2/2,0.01);
        add=add/2;
    end
    if (tk1~=Le) && (length(Code1)<1000) 
        tk1=min(tk+4000,Le);
        Code1=[Code1,Code(tk+1:tk1)];
        tk=tk1; 
    end
end




