function [Lw,Up,cnew]=enc(value,Lw1,Up1,pb)
% encoding for arithmetic coding
pcum=[0,cumsum(pb)];
pcum=pcum/pcum(end);
p1=pcum(value);
p2=pcum(value+1);
lu=Up1-Lw1;
Lw=Lw1+lu*p1; 
Up=Lw1+lu*p2;
cnew=[];
while (Lw>=0.5) || (Up<=0.5)
    if (Lw>=0.5)
        Lw=Lw*2-1;   Up=Up*2-1;
        cnew=[cnew,'1'];
    else
        Lw=Lw*2;     Up=Up*2;
        cnew=[cnew,'0'];
    end
end
    


