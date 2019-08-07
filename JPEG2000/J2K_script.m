%#################################### README ###################################################################
% This code is an implementation of encoding & decoding JPEG2000 (CX,D) pairs using the proposed method
% To run this Code:
%       1. Make sure you are in the right directory with the according files.
%       2. Change the directory of your pre-generated (CX,D) pairs from JPEG2000 coding passes.
%		   (CX,D) pairs should be of the size : <2xN double>, CX in the 1st row, and D in the 2nd row.
%		   some examples are included in the directory: mat.
% --------------------------------------------------------------------------------------------------------------
% Last update: 05/23/2017
% Contact: b01901091@ntu.edu.tw
%###############################################################################################################
clear;clc;
addpath .\mat;
load lenaYCXD.mat;  YCXD=CXD;
load lenaCbCXD.mat; CbCXD=CXD;
load lenaCrCXD.mat; CrCXD=CXD;
% Encode & decode using the proposed method
Ycode=J2K_enco(YCXD);
Cbcode=J2K_enco(CbCXD);
Crcode=J2K_enco(CrCXD);
lena_b = length(Ycode) + length(Cbcode) + length(Crcode);
fprintf('code length = %g bits, %g bpp (originally bpp = 24 bits) \n', lena_b, lena_b/512/512);
YCXD_new=J2K_deco(Ycode,YCXD(1,:));
CbCXD_new=J2K_deco(Cbcode,CbCXD(1,:));
CrCXD_new=J2K_deco(Crcode,CrCXD(1,:));
if isequal(YCXD,YCXD_new) && isequal(CbCXD,CbCXD_new) && isequal(CrCXD,CrCXD_new)
    fprintf('valid\n');
else
    fprintf('invalid\n');
end