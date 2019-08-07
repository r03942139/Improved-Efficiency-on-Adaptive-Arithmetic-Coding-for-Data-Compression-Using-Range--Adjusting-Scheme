function [vn,Code,Lw,Up,Lw2,Up2]=dcd(Code1,Lw1,Up1,Lw2,Up2,pb)
% decoding for arithmetic coding
% Code1=C1; Lw2=10/16; Up2=11/16; pb=pm; 
% Lw1=Lw; Up1=Up; 
pcum=[0,cumsum(pb)];
pcum=pcum/pcum(end);
lu=Up1-Lw1;
p1=Lw1+lu*pcum;
Ln=length(pb);
bt=1;
fd=0; fst=0;
while~fd
    vn=fst+find(p1(fst+1:Ln+1)>Lw2,1,'first')-1;    
    fst=vn-1;    
    if (vn<=Ln) && (p1(vn+1)>=Up2)
        fd=1;
    elseif Code1(bt)=='1';           
        Lw2=Lw2+(Up2-Lw2)/2;   bt=bt+1;
    else
        Up2=Lw2+(Up2-Lw2)/2;   bt=bt+1;
    end    
end
Lw=p1(vn);    Up=p1(vn+1); 
Code=Code1(bt:end);
while (Lw>=0.5) || (Up<=0.5)
    if (Lw>=0.5)
        Lw=Lw*2-1;   Up=Up*2-1;
        Lw2=Lw2*2-1; Up2=Up2*2-1;
    else
        Lw=Lw*2;     Up=Up*2;
        Lw2=Lw2*2;   Up2=Up2*2;
    end
end
    


