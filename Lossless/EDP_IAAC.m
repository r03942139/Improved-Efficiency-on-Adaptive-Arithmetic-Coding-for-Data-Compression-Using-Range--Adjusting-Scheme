function [Code,bpp]=EDP_IAAC(A1)
%%  Edge-Directed Prediction  
%% Standard Code of the Proposed IAAC Algorithm for EDP Lossless Encoding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
A1=double(A1);
%% Initialization
[H,W]=size(A1);
MN = H*W;   
img_src = reshape( A1',1, [] );            
E= zeros( H,W);   % prediction residue image
Ie=zeros(H,W);       % predicted image
%% Training window       
% Neighbors
N=10;
Ny=[ 0, 0, -1, -1, -1, -1, -1, -2, -2, -2 ];
Nx=[ -1,-2, -2, -1,  0, 1, 2, -1,  0, 1 ];
xs1=-min(Nx); xs2=max(Nx); ys=-min(Ny);
T=min(N,7);
TT= 2*T*(T+1);
Tx=[kron(ones(1,T),[-T:T]),[-T:-1]];
Ty=[kron([-T:-1],ones(1,2*T+1)),zeros(1,T)];
Txy = Ty*W + Tx;
th=2;
Tn1=T+ys;
T2w=kron([1:T]*W,ones(1,T+1))-kron(ones(1,T),[0:T]);
ts=length(T2w);
T1w=[T2w'+W,T2w'+2*W];
T3w=[W,W*2];
WD1=W-xs2;
WT1=W-T+1-xs2;
T2=T+ys;
T3=T+xs1;
Tm=(Ty'*ones(1,N)+ones(TT,1)*Ny)*W+Tx'*ones(1,N)+ones(TT,1)*Nx;
Nxy=Ny*W+Nx;
Mx=[-1,0,-1]; My=[-1,-1,0]; Mxy=Mx+My*W;
Nm=length(Mx);
Tm2=(Ty'*ones(1,Nm)+ones(TT,1)*My)*W+Tx'*ones(1,Nm)+ones(TT,1)*Mx;
flag1=1;  em1=0;
flag3=ones(1,xs1+W-WD1+1);
em3=zeros(1,xs1+W-WD1+1);  L3=zeros(1,xs1+W-WD1+1);
x12=[2:xs1+1,WD1:W];
%% Initialization for Parameters for AAC
bpp=zeros(1,5);
s1=0.7;
sca=s1;
ac=0.8;
b=2.4;  
cc=14;   
cc1=9;
cc2=cc1*2;
ck=cc*cc1;
map=[1,2,3,3,4,4,5,6];
th0=[0.5,2.6*1.75.^[1:8]];      
ln2=0.2.^([cc1-.5:-1:.5]/cc1).*0.06.^([.5:cc1-.5]/cc1)*0.6;
dx1=[1:cc]/cc*8;
dx2=floor(dx1)+1;
dx3=dx1-dx2+1;
th1=th0(dx2(1:cc-1)).^(1-dx3(1:cc-1)).*th0(dx2(1:cc-1)+1).^dx3(1:cc-1);
mz1=0.6+1.3.^[1:cc]*0.05;
sz1=2+2*1.4.^[1:cc]*0.11;
mz=kron(mz1,[-(cc1-1)/2:(cc1-1)/2]*3.3/cc1);
sz=kron(sz1,ones(1,cc1));
xx=[-255:255];
eps=zeros(256,511);
eps(1,256)=1;
for e21=1:255
    sd=abs(e21)/8;
    eps(e21+1,:)=exp(-abs(xx-e21)/sd);
    eps(e21+1,:)=eps(e21+1,:)/sum(eps(e21+1,:));
end
eps=[eps(256:-1:2,511:-1:1);eps];
Tb1=zeros(ck,511);
for b1=1:ck
    Tb1(b1,:)=exp(-(xx-mz(b1)).^2/2/sz(b1)^2);
    Tb1(b1,:)=Tb1(b1,:)/sum(Tb1(b1,:))*10000;
    Tb1(b1,:)=max(Tb1(b1,:),0.1);
end
cz1=kron(ones(1,cc),[0:cc1-1]);
cz=ones(cc1-1,1)*cc1*floor([0:ck-1]/cc1)+1;
for b2=1:cc1-1;
    cz2=mod(cz1+b2,cc1);
    cz(b2,:)=cz(b2,:)+cz2;
end
mz2=round(mz);
xn2=ones(cc1-1,1)*[2:512];
AD=500*ones(1,ck);  % adjustable  2
Up=1;
Lw=0;
%%  Boundary Prediction by DPCM
mul1=[3,1.5,ones(1,W-3),1.5];
mul2=[1.5,ones(1,W-2),1.5];
C1=zeros(H,W);
Ctx=zeros(H,W);
T1=Tb1;
Ic=zeros(H,W);
E(1,1)= img_src(1);               % 1st row
Code=[];
Code2=dec2bin(img_src(1),8);           % Code
% L1=zeros(H,W);  U1=zeros(H,W);  
% ep11=zeros(H,W); ep11(1,1)=0;
for m1=1:H
    if m1==1
        dh=[0,0,abs(A1(1,2:W-1)-A1(1,1:W-2))]*3;
        dv=zeros(1,W);     
        dh=dh*sca;        
        Ic(1,1)=128;
        Ic(1,2:W)=A1(1,1:W-1);   
        Ie(1,2:W)=A1(1,1:W-1);
        E(1,2:W)= img_src( 2:W) - Ie(1,2:W); 
    else
        % the first element
        Ie1=zeros(1,W);
        es=zeros(1,W);
        if m1<=Tn1                    
            Ie1(1)=A1(m1-1,1);  es(1)=A1(m1,1)-Ie1(1);
        else
            x=(m1-1)*W+1;
            if flag1
                C=img_src(x-T1w);
                A= C' * C;
                if min(abs(eig(A)))>0.00001
                    B= C' *img_src(x-T2w)';
                    ac1= A \ B;
                    em1=1;
                end
            end
            if em1
                Ie1(1)=  round(img_src( x-T3w)*ac1);
            else
                Ie1(1)=  round(img_src( x-W));
            end
            if(Ie1(1) <0 || Ie1(1)>255 )  Ie1(1)= round(img_src( x-W));  end
            es(1)=img_src( x)-Ie1(1);
            flag1=(abs(es(1))>=th)||(~em1);
        end
        % the 2nd near-border part
        if m1 <=T+1
            if m1<=ys+1
                xc=[2:T+1,W-T+1:W];
            else
                xc=[2:xs1+1,WD1:W];
            end
            xa=(m1-1)*W+xc;
            pn1  = img_src(xa-W);  % north      
            pw1  = img_src(xa-1);  % west
            pnw1 = img_src(xa-W-1); 
            min1= min(pn1,pw1);      max1= max(pn1,pw1);
            cs=(pnw1>=max1)+((pnw1<=min1)&(pnw1<max1))*2;
            Ie1(xc)=round(min1.*(cs==1)+max1.*(cs==2)+(pn1+pw1-pnw1).*(cs==0));
            es( xc)= img_src(xa)-Ie1(xc);   
        end
        % 2nd and 3rd row
        if m1 <=ys+1
            f1=find(Ty>=2-m1);
            Tm3=Tm2(f1,:);
            Txy3=Txy(f1);
            flag2=1;   em2=0;
            for x11=T+2:W-T
                x=(m1-1)*W+x11;
                if flag2
                    C=img_src(x+Tm3);            
                    A= C' * C;
                    if min(abs(eig(A)))>0.00001
                        B= C' *img_src(x+Txy3)';
                        ac2= A \ B;
                        em2=1;
                    end
                end
                if em2
                    pred=  img_src( x+Mxy)*ac2;
                else
                    xw=img_src( x-1);   xn=img_src( x-W);    xwn=img_src( x-W-1);
                    max1=max(xw,xn);    min1=min(xw,xn);
                    if xwn>=max1
                        pred=min1;
                    elseif xwn<=min1
                        pred=max1;
                    else
                        pred=xw+xn-xwn;
                    end
                end
                if( pred <0 || pred>255 )  pred= mean(img_src(x+Mxy));  end
                Ie1(x11)=round(pred);
                es(x11)=img_src( x)-Ie1(x11);
                flag2=(abs(es(x11))>=th)||(~em2);
            end   
        end
        % 2nd, 3rd, and the last three columns
        if m1>=T+2
            for a3=1:xs1+W-WD1+1
                x13=x12(a3);
                f1=find((Tx>=2-x13)&(Tx<=W-x13));
                Tm3=Tm2(f1,:);
                Txy3=Txy(f1);
                x=(m1-1)*W+x13;
                if flag3(a3)
                    C=img_src(x+Tm3);            
                    A= C' * C;
                    if m1==T+2
                        [sa1,sa2]=size(A);
                        L3(a3)=sa2;
                    end
                    if min(abs(eig(A)))>0.00001
                        B= C' *img_src(x+Txy3)';
                        ac3(1:L3(a3),a3)= A \ B;
                        em3(a3)=1;
                    end
                end
                if em3(a3)
                    pred=  img_src( x+Mxy)*ac3(1:L3(a3),a3);
                else
                    xw=img_src( x-1);   xn=img_src( x-W);    xwn=img_src( x-W-1);
                    max1=max(xw,xn);    min1=min(xw,xn);
                    if xwn>=max1
                        pred=min1;
                    elseif xwn<=min1
                        pred=max1;
                    else
                        pred=xw+xn-xwn;
                    end
                end
                if( pred <0 || pred>255 )  pred= mean(img_src(x+Mxy));  end
                Ie1(x13)=round(pred);
                es(x13)=img_src( x)-Ie1(x13);
                flag3(a3)=(abs(es(x13))>=th)||(~em3(a3));
            end
        end           
              %% Least-square Prediction (LSP)
        if m1>ys+1    
            xe=[];
            flag=1;  em=0; 
            for n1=xs1+2:WD1-1   
                ds=0;         
                x=(m1-1)*W+n1;
                if( m1<=T2  ||  n1<=T3  ||  n1>=WT1 )  % another special case
                    if flag
                        yz=min(m1-ys-1,T);
                        xz1=min(n1-xs1-1,T);
                        xz2=min(WD1-n1,T);
                        f1=find((Ty>=-yz)&(Tx>=-xz1)&(Tx<=xz2));
                        Txy1=Txy(f1);
                        Tm1=Tm(f1,:);
                        % 2nd-sub border
                        C=img_src(x+Tm1);
                        A= C' * C;     
                        if min(abs(eig(A)))>0.00001
                            B= C' * img_src(x+Txy1)';
                            a= A \ B;    
                            em=1;
                        end
                    end
                else       
                    if flag
                        C=img_src(x+Tm);
                        A= C' * C;
                        if min(abs(eig(A)))>0.00001
                            B= C' *img_src(x+Txy)';
                            a= A \ B;
                            em=1;
                        end
                    end
                end
                if em
                    pred=  img_src( x+Nxy)*a;
                    if( pred <0 || pred>255 )  pred= mean(img_src(x+Nxy));  end
                    Ie1(n1)=round(pred);
                    es(n1)=img_src( x)-Ie1(n1);
                    flag=(abs(es(n1))>=th);
                % x/MN                
                else
                    xe=[xe,x];
                    flag=1;
                end    
            end
            xw=img_src( xe-1);    xn=img_src( xe-W);     xwn=img_src( xe-W-1);
            max1=max(xw,xn);      min1=min(xw,xn);
            cs=(xwn>=max1)+((xwn<=min1)&(xwn<max1))*2;
            pred=min1.*(cs==1)+max1.*(cs==2)+(xn+xw-xwn).*(cs==0);
            xe1=xe-(m1-1)*W;
            Ie1(xe1)=round(pred);
            es(xe1)= img_src(xe)-Ie1(xe1);  
        end        
        Ie(m1,:)=Ie1;
        E(m1,:)=es;
        dh1=[0,abs(A1(m1-1,2:W)-A1(m1-1,1:W-1))];
        dh2=[dh1(2:W),0];
        dh3=[0,0,abs(A1(m1,2:W-1)-A1(m1,1:W-2))];
        dh=(dh1+dh2+dh3).*mul1;
        if m1==2
            dv=[0,abs(A1(2,1:W-1)-A1(1,1:W-1))]*3;
        else
            dv1=abs(A1(m1-1,:)-A1(m1-2,:));
            dv2=[dv1(2:W),0];
            dv3=[0,abs(A1(m1,1:W-1)-A1(m1-1,1:W-1))];
            dv=(dv1+dv2+dv3).*mul2;
        end
        dh=dh*sca;
        dv=dv*sca;
        I4=(A1(m1,1:W-1)+A1(m1-1,2:W))/2+(A1(m1-1,[3:W,W])-A1(m1-1,[1:W-2,W]))*0.15;
        dvh=dv(2:W)-dh(2:W);
        Cas=ones(1,W-1)+(dvh>=-80)+(dvh>=-32)+(dvh>=-8)+(dvh>8)+(dvh>32)+(dvh>80);
        In=A1(m1-1,2:W);  
        Iw=A1(m1,1:W-1);
        A11=(Cas==1).*In+(Cas==7).*Iw+(Cas==4).*I4+(Cas==2).*(I4+In)/2+...
            (Cas==3).*(I4*3+In)/4+(Cas==5).*(I4*3+Iw)/4+(Cas==6).*(I4+Iw)/2;
        Ic(m1,1)=A1(m1-1,1);
        Ic(m1,2:W)=min(max(round(A11),0),255);
    end
%    Step 3: Context Modeling for Gradients
    ew1=E(m1,:);
    ew2=abs(ew1);
    if m1==1
        ew=[0,0,ew2(2:W-1)];
    else
        ew=[abs(E(m1-1,1)),(abs(E(m1-1,2:W))+ew2(1:W-1))/2];
        if m1==2
            ew(1)=0;
        end
    end
    delta=ac*(dh+dv)+b*ew;
    if m1==1
        delta(3:W)=delta(3:W)+dh(3:W);
    elseif m1==2
        delta(1)=delta(1)+dh(1);
    end
    Ce1=ones(1,W);
    for b1=1:cc-1
        Ce1=Ce1+(delta>th1(b1));
    end
    C1(m1,:)=Ce1;
    % Step 4: Texture Context
    Ice=Ic(m1,:)-Ie(m1,:);
    C2=1+(Ice>=-7)+(Ice>=-3)+(Ice>=-1)+(Ice>=0)+(Ice>0)+(Ice>1)+(Ice>3)+(Ice>7);
    Ctx(m1,:)=(Ce1-1)*cc1+C2;
    % Step 5: Adaptive Arithmetic Coding
    n0=1+(m1==1);  
    for n=n0:W
        C=Ctx(m1,n);
        T12=T1(C,:);
        T13=T12;      % local frequency table
        if n>1
            kt=256+A1(m1,n-1)-Ie(m1,n);
            T13(kt)=T13(kt)*1.05; 
        end
        if m1>1
            kt=256+A1(m1-1,n)-Ie(m1,n);
            T13(kt)=T13(kt)*1.05; 
        end
        ts1=255-Ie(m1,n);
        T13=T13(ts1+1:ts1+256); 
        E2=E(m1,n);
        % entropy=entropy+log2(Ts/T13(256+E2));          
        [Lw,Up,cnew]=enc(256+E2-ts1,Lw,Up,T13);
        % L1(m1,n)=Lw;  U1(m1,n)=Up; 
        Code2=[Code2,cnew];
        % ep11(m1,n)=length(Code)+length(Code2)-log2(Up-Lw);
        if E2==0
            ep1=zeros(1,511);
            ep1(256)=1;
        else
            ep1=eps(256+E2,:);
        end
        T12=T12+ep1*AD(C);
        if sum(T13)+AD(C)>1000000
              T12=max(T12/2,0.1);
              AD(C)=AD(C)*0.5;
        end
        T1(C,:)=T12;  % mutual learning
        ep31=ep1*0.6;
        if C-cc1>0
            T1(C-cc1,:)=T1(C-cc1,:)+ep31*AD(C-cc1);
        end
        if C+cc1<=ck
            T1(C+cc1,:)=T1(C+cc1,:)+ep31*AD(C+cc1);
        end
        ofn=mz2(cz(:,C))'-mz2(C);
        cd1=abs(C-cz(:,C));
        ep2=[0,ep1,0];
        xn1=min(max(xn2-ofn*ones(1,511),1),513);
        T1(cz(:,C),:)=T1(cz(:,C),:)+diag(ln2(cd1).*AD(cz(:,C)))*ep2(xn1);       
%         if (m1==2) && (n==1)
%             [Ie(2,1),C,delta(1),dh(1),dv(1),ew(1)]
%             return
%         end
    end 
    AD=AD*1.02;
    sgm=max(mean((ew1-mean(ew1)).^2)^.5,1);    
    sca=s1*1.55^(log2(sgm)-2.8);
    if length(Code2)>2000
        Code=[Code,Code2];
        Code2=[];
    end
    row=m1
end
Code=[Code,Code2];
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
E1=E;
bpp=length(Code)/MN;
% save('LCode_Lena.mat','Code');
% save('LCode_Lena_test.mat','L1','U1');



