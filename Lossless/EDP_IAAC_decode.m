function A2=EDP_IAAC_decode(Code,H,W)
%%  Edge-Directed Prediction  
%% Standard Code of the Proposed IAAC Algorithm for EDP Lossless Decoding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%% Initialization
% clear 
% load('LCode_Lena.mat');
L=length(Code);
H=512; W=512; 
A2=zeros(W,H);
MN = H*W;      
E= zeros(W,H);      % prediction residue image
Ie=zeros(W,H);       % predicted image
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
Lw=0;  Up=1;
Lw2=0; Up2=1;
%%  Boundary Prediction by DPCM
mul1=[3,1.5,ones(1,W-3),1.5];
mul2=[1.5,ones(1,W-2),1.5];
C1=zeros(W,H);
Ctx=zeros(W,H);
T1=Tb1;
Ic=zeros(W,H);  Ic(1,1)=128;
ckn=min(4000,L)+8;
Code5=Code(9:ckn);
A2(1,1)=bin2dec(Code(1:8));           % Case 1 
nc51=zeros(1,W);   nc51([2:(T+1),(W-T+1):W])=1;
nc52=zeros(1,W);   nc52([2:(xs1+1),WD1:W])=1;
nc6=zeros(1,W);    nc6([(T+2):(W-T)])=1;
nc7=zeros(1,W);    nc7(x12)=[1:length(x12)];
nc8=zeros(1,W);    nc8([(xs1+2):T3,WT1:(WD1-1)])=1;
nc9=zeros(1,W);    nc9((xs1+2):(WD1-1))=1;
% L2=zeros(H,W);     U2=zeros(H,W);
for m1=1:H
    n0=1+(m1==1);  
    sca3=3*sca;
    flag=1;  em=0; 
    ek=0;
    for n=n0:W
        x=(m1-1)*W+n;
        if m1==1                       % Case 2
            if n==2
                dh=0;
            else
                dh=abs(A2(n-1,1)-A2(n-2,1))*sca3;
            end
            dv=0;   
            Ic(n,1)=A2(n-1,1);   
            Ie(n,1)=A2(n-1,1);
        else
            % the first element
            % Ie1=zeros(1,W);  
            if n==1       
                if m1<=Tn1            % Case 3      
                    Ie(1,m1)=A2(1,m1-1);                  
                else                  % Case 4
                    if flag1
                        C=A2(x-T1w);
                        A= C' * C;
                        if min(abs(eig(A)))>0.00001
                            B= C' *A2(x-T2w)';
                            ac1= A \ B;
                            em1=1;
                        end
                    end
                    if em1
                        Ie(1,m1)=  round(A2( x-T3w)*ac1);
                    else
                        Ie(1,m1)=  round(A2(1,m1-1));
                    end
                    if(Ie(1,m1) <0 || Ie(1,m1)>255 )  Ie(1,m1)= round(A2(1,m1-1));  end                    
                end
            end
            % the 2nd near-border part
            if ((m1 <=ys+1 ) && nc51(n))||((m1 <=T+1 )&& nc52(n))   % Case 5  
                pn1  = A2(x-W);  % north      
                pw1  = A2(x-1);  % west
                pnw1 = A2(x-W-1); 
                min1= min(pn1,pw1);      max1= max(pn1,pw1);
                cs=(pnw1>=max1)+((pnw1<=min1)&(pnw1<max1))*2;
                Ie(n,m1)=round(min1*(cs==1)+max1*(cs==2)+(pn1+pw1-pnw1)*(cs==0));
            end
            % 2nd and 3rd row
            if m1 <=ys+1                         % Case 6
                if n==T+2
                    f11=find(Ty>=2-m1);
                    Tm31=Tm2(f11,:);
                    Txy31=Txy(f11);
                    flag2=1;   em2=0;
                end
                if nc6(n)
                    if flag2
                        C=A2(x+Tm31);            
                        A= C' * C;
                        if min(abs(eig(A)))>0.00001
                            B= C' *A2(x+Txy31)';
                            ac2= A \ B;
                            em2=1;
                        end
                    end
                    if em2
                        pred=  A2( x+Mxy)*ac2;
                    else
                        xw=A2( x-1);   xn=A2( x-W);    xwn=A2( x-W-1);
                        max1=max(xw,xn);    min1=min(xw,xn);
                        if xwn>=max1
                            pred=min1;
                        elseif xwn<=min1
                            pred=max1;
                        else
                            pred=xw+xn-xwn;
                        end
                    end
                    if( pred <0 || pred>255 )  pred= mean(A2(x+Mxy));  end
                    Ie(n,m1)=round(pred);
                end   
            end
            % 2nd, 3rd, and the last three columns
            if (m1>=T+2) && nc7(n)                               % Case 7
                a3=nc7(n);
                f12=find((Tx>=2-n)&(Tx<=W-n));
                Tm32=Tm2(f12,:);
                Txy32=Txy(f12);
                if flag3(a3)
                    C=A2(x+Tm32);            
                    A= C' * C;
                    if m1==T+2
                        [sa1,sa2]=size(A);
                        L3(a3)=sa2;
                    end
                    if min(abs(eig(A)))>0.00001
                        B= C' *A2(x+Txy32)';
                        ac3(1:L3(a3),a3)= A \ B;
                        em3(a3)=1;
                    end
                end
                if em3(a3)
                    pred=  A2( x+Mxy)*ac3(1:L3(a3),a3);
                else
                    xw=A2( x-1);   xn=A2( x-W);    xwn=A2( x-W-1);
                    max1=max(xw,xn);    min1=min(xw,xn);
                    if xwn>=max1
                        pred=min1;
                    elseif xwn<=min1
                        pred=max1;
                    else
                        pred=xw+xn-xwn;
                    end
                end
                if( pred <0 || pred>255 )  pred= mean(A2(x+Mxy));  end
                Ie(n,m1)=round(pred);                
            end           
                    %% Least-square Prediction (LSP)
            if (m1>ys+1) && nc9(n)                            % Cases 8, 9                  
                ds=0;  
                if( m1<=T2  ||  nc8(n)  )                   % Case 8
                    if flag
                        yz=min(m1-ys-1,T);
                        xz1=min(n-xs1-1,T);
                        xz2=min(WD1-n,T);
                        f1=find((Ty>=-yz)&(Tx>=-xz1)&(Tx<=xz2));
                        Txy1=Txy(f1);
                        Tm1=Tm(f1,:);
                        C=A2(x+Tm1);
                        A= C' * C;     
                        if min(abs(eig(A)))>0.00001
                            B= C' * A2(x+Txy1)';
                            a= A \ B;    
                            em=1;
                        end
                    end
                else                                       % Case 9
                    if flag
                        C=A2(x+Tm);
                        A= C' * C;
                        if min(abs(eig(A)))>0.00001
                            B= C' *A2(x+Txy)';
                            a= A \ B;
                            em=1;
                        end
                    end
                end
                if em
                    pred=  A2( x+Nxy)*a;
                    if( pred <0 || pred>255 )  pred= mean(A2(x+Nxy));  end
                    Ie(n,m1)=round(pred);
                    ek=1;
                else
                    xw=A2( x-1);    xn=A2( x-W);     xwn=A2( x-W-1);
                    max1=max(xw,xn);      min1=min(xw,xn);
                    cs=(xwn>=max1)+((xwn<=min1)&(xwn<max1))*2;
                    pred=min1*(cs==1)+max1*(cs==2)+(xn+xw-xwn)*(cs==0);
                    Ie(n,m1)=round(pred);
                    flag=1;                  
                end               
            end 
            if n==1
                dh1=0;
            else
                dh1=abs(A2(n,m1-1)-A2(n-1,m1-1));
            end
            if n < W
                dh2=abs(A2(n+1,m1-1)-A2(n,m1-1));
            else
                dh2=0;
            end
            if n<=2
                dh3=0;
            else
                dh3=abs(A2(n-1,m1)-A2(n-2,m1));
            end
            dh=(dh1+dh2+dh3)*mul1(n);
            if m1==2
                if n==1
                    dv=0;
                else
                    dv=abs(A2(n-1,2)-A2(n-1,1))*3;
                end
            else
                dv1=abs(A2(n,m1-1)-A2(n,m1-2));
                if n<W
                    dv2=abs(A2(n+1,m1-1)-A2(n+1,m1-2));
                else
                    dv2=0;
                end
                if n==1
                    dv3=0;
                else
                    dv3=abs(A2(n-1,m1)-A2(n-1,m1-1));
                end
                dv=(dv1+dv2+dv3)*mul2(n);
            end
            dh=dh*sca;
            dv=dv*sca;    
            % Step 2
            if n==1
                Ic(1,m1)=A2(1,m1-1);
            else
                if n ~=W
                    I4=(A2(n-1,m1)+A2(n,m1-1))/2+(A2(n+1,m1-1)-A2(n-1,m1-1))*0.15;
                else
                    I4=(A2(W-1,m1)+A2(W,m1-1))/2;
                end
                dvh=dv-dh;
                Cas=1+(dvh>=-80)+(dvh>=-32)+(dvh>=-8)+(dvh>8)+(dvh>32)+(dvh>80);
                In=A2(n,m1-1);  
                Iw=A2(n-1,m1);
                A11=(Cas==1)*In+(Cas==7)*Iw+(Cas==4)*I4+(Cas==2)*(I4+In)/2+...
                    (Cas==3)*(I4*3+In)/4+(Cas==5)*(I4*3+Iw)/4+(Cas==6)*(I4+Iw)/2;
                Ic(n,m1)=min(max(round(A11),0),255);
            end        
        end      
       %  Step 3: Context Modeling for Gradients
        % ew1=E(m1,:);       % ew2=abs(ew1);
        if m1==1
            if n==2
                ew=0;
            else
                ew=abs(E(n-1,1));
            end
        else
            if n==1
                if m1==2
                    ew=0;
                else
                    ew=abs(E(1,m1-1));
                end
            else
                ew=(abs(E(n,m1-1))+abs(E(n-1,m1)))/2;            
            end
        end
        delta=ac*(dh+dv)+b*ew;
        if ((m1==1) && (n>2)) || ((m1==2) && (n==1))
            delta=delta+dh;
        end
        Ce1=1+sum(delta>th1);        
        C1(n,m1)=Ce1;
        % Step 4: Texture Context
        Ice=Ic(n,m1)-Ie(n,m1);
        C2=1+(Ice>=-7)+(Ice>=-3)+(Ice>=-1)+(Ice>=0)+(Ice>0)+(Ice>1)+(Ice>3)+(Ice>7);
        Ctx(n,m1)=(Ce1-1)*cc1+C2;
        % Step 5: Adaptive Arithmetic Coding
        C=Ctx(n,m1);
        T12=T1(C,:);
        T13=T12;      % local frequency table
        if n>1
            kt=256+A2(n-1,m1)-Ie(n,m1);
            T13(kt)=T13(kt)*1.05; 
        end
        if m1>1
            kt=256+A2(n,m1-1)-Ie(n,m1);
            T13(kt)=T13(kt)*1.05; 
        end
        ts1=255-Ie(n,m1);
        T13=T13(ts1+1:ts1+256); 
        [vn,Code5,Lw,Up,Lw2,Up2]=dcd(Code5,Lw,Up,Lw2,Up2,T13);
        %  L2(m1,n)=Lw;   U2(m1,n)=Up;  
        E2=vn-256+ts1;
        E(n,m1)=E2;
        A2(n,m1)=Ie(n,m1)+E2;
        % [Lw,Up,cnew]=enc(256+E2-ts1,Lw,Up,T13);
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
        if (n==1)&&(m1>Tn1)
            flag1=(abs(E2)>=th)||(~em1);
        elseif (m1 <=ys+1) && (m1>1) && nc6(n) 
            flag2=(abs(E2)>=th)||(~em2);  
        elseif (m1>=T+2) && nc7(n)
            flag3(a3)=(abs(E2)>=th)||(~em3(a3));
        elseif ek
            flag=(abs(E2)>=th);
        end
        if length(Code5)<1000
            ckn1=min(ckn+4000,L);
            Code5=[Code5,Code(ckn+1:ckn1)];
            ckn=ckn1;
        end
%         if (abs(L2(m1,n)-L1(m1,n))>0.0001) && (abs(U2(m1,n)-U1(m1,n))>0.0001)
%             [m1,n]
%             [Ie(n,m1),C,delta,dh,dv,ew]
%             return
%         end
    end
    AD=AD*1.02;
    ew1=E(:,m1);
    if m1==1  ew1(1)=A2(1,1);  end
    sgm=max(mean((ew1-mean(ew1)).^2)^.5,1);    
    sca=s1*1.55^(log2(sgm)-2.8);
    row=m1
end
A2=A2.';
% if sum(sum(abs(A2-A1)))==0
%     'Successful!'
% end


