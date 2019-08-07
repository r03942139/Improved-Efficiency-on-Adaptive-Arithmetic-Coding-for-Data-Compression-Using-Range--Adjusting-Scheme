A=imread('C:\Users\DJJ\Documents\Pic\Lena_gray_512.bmp');
A=double(A);
[M,N]=size(A);
M1=M/8;
N1=N/8;
C=cos(pi*[0:7]'*[.5:7.5]/8)/2;
C(1,:)=C(1,:)/2^.5;
C1=C';
dc=zeros(1,M1*N1);
kn=0;
for a1=1:M1
    vn=zeros(1,N1);
    x1=a1*8+[-7:0];
    for a2=1:N1
        A1=A(x1,a2*8+[-7:0]);
        cf=C*A1*C1;
        vn(a2)=round(cf(1,1)/16);
    end
    dc(kn+[1:N1])=vn;
    kn=kn+N1;
end

Code = DC_IAAC(dc,M,N);
length(Code)

dc2 = DC_IAAC_decode(Code,M,N);

if sum(abs(dc2-dc))==0
    'Successful!'
end

