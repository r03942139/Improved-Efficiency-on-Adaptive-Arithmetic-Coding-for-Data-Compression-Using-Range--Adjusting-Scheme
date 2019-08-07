function Code = DC_IAAC(dc,M,N)
% Standard code using the proposed IAAC algorithm to encode the DC terms 
dcd=diff(dc);
M1=M/8;  N1=N/8;  MN=M1*N1;
Lw=0; Up=1;
[Lw,Up,Code]=enc(dc(1)+1,Lw,Up,ones(1,129)/129);
tb2=exp(-0.2*[0:128]);
tb2=tb2/sum(tb2)*1000;
tb2=max(tb2,1);
dc1=dc(1);
add=80;
for a7=1:(MN-1)
    d2=abs(dcd(a7));
    d1=d2+1;
    dm=max(128-dc1,dc1);
    dm1=min(128-dc1,dc1);
    tb22=tb2(1:dm+1);
    if mod(a7,N)==0
        tb22=tb22.*(1+log2(1:dm+1)*1);
    end
    tb22(dm1+2:dm+1)=tb22(dm1+2:dm+1)*.8;   
    [Lw,Up,cnew]=enc(d1,Lw,Up,tb22);
    Code=[Code,cnew];
    if (d1~=1)&&(dc1-d2>=0)&&(dc1+d2<=128)
        [Lw,Up,cnew]=enc(1+(sign(dcd(a7))==1),Lw,Up,[.5,.5]);
        Code=[Code,cnew];
    end
    if d2==0
        ad1=[add,zeros(1,128)];
    else
        d3=abs(([1:128]+1)/(d2+1)-1);
        ad1=[0,exp(-7*d3.^2)];
        ad1=ad1/sum(ad1)*add;        
    end
    tb2=tb2+ad1;
    dc1=dc1+dcd(a7);  
    add=add*1.02;            % adjusting step
    if sum(tb2)>100000
        tb2=max(tb2/2,0.01);
        add=add/2;
    end
end
fn=(Up<1)||(Lw>0);
bn=0;
while fn
    bn=bn+1;
    Lw=Lw*2;
    Up=Up*2;    
    fn=Up<ceil(Lw)+1;
end
if bn
    Code=[Code,dec2bin(ceil(Lw),bn)];
end

