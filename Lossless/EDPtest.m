A1=double(imread('C:\Users\DJJ\Documents\Pic\gray512\Baboon.bmp'));
[H,W]=size(A1);
[Code,bpp]=EDP_IAAC(A1);
A2=EDP_IAAC_decode(Code,H,W);
if sum(sum(abs(A2-A1)))==0
    'Successful!'
end
bpp
