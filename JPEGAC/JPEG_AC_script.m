%#################################### README ###################################################################
% This code is an implementation of encoding & decoding JPEG AC coefficients using the proposed method
% To run this Code:
%       1. Make sure you are in the right directory with the according files.
%       2. Execute the script
% --------------------------------------------------------------------------------------------------------------
% Last update: 05/23/2017
% Contact: b01901091@ntu.edu.tw
%###############################################################################################################


clear;clc;
tic
result = zeros(8,1);
for img = 1:8
if img==1
    A=imread([pwd,'\img\lena.bmp']);
elseif img==2
    A=imread([pwd,'\img\baboon.bmp']);
elseif img==3
    A=imread([pwd,'\img\peppers.bmp']);
elseif img==4
    A=imread([pwd,'\img\airplane.bmp']);
elseif img==5
    A=imread([pwd,'\img\tiffany.bmp']);
elseif img==6
    A=imread([pwd,'\img\splash.jpg']);
elseif img==7
    A=imread([pwd,'\img\sailboat.jpg']);
elseif img==8
    A=imread([pwd,'\img\house1.bmp']);
end
A=double(A);
%A = floor(bitsra(A,1));
load JPEGtable.mat
%coding
[M,N,c3]=size(A);
Y=0.299*A(:,:,1)+0.587*A(:,:,2)+0.114*A(:,:,3);
Cb=-0.169*A(:,:,1)-0.331*A(:,:,2)+0.5*A(:,:,3);
Cr=0.5*A(:,:,1)-0.419*A(:,:,2)-0.081*A(:,:,3);
Cb1=Cb(1:2:M,1:2:N)+127.5;
Cr1=Cr(1:2:M,1:2:N)+127.5;
M1=M/8;   N1=N/8;
M2=M1/2;  N2=N1/2;
A1=Y; A11=Cb1;  A12=Cr1;
x1=[-7:0];  y1=[-7:0];
for a1=1:M1
    x1=x1+8;
    A1(x1,:)=C*A1(x1,:);
    if a1<=M2
        A11(x1,:)=C*A11(x1,:);
        A12(x1,:)=C*A12(x1,:);
    end
end
for a2=1:N1
    y1=y1+8;
    A1(:,y1)=A1(:,y1)*C1;
    if a2<=N2
        A11(:,y1)=A11(:,y1)*C1;
        A12(:,y1)=A12(:,y1)*C1;
    end
end
Qc=1./Q; Q1c=1./Q1;
YAC = zeros(M1,N1,63);  CbAC = zeros(M2,N2,63); CrAC = zeros(M2,N2,63);
x1=(-7:0);  
for a3=1:M1
    x1=x1+8;
    y1=(-7:0); 
    for a4=1:N1
        y1=y1+8;
        A2=round(A1(x1,y1).*Qc);
        YAC(a3,a4,:) = zigzag(A2);
        if (a3<=M2)&&(a4<=N2)
            A21=round(A11(x1,y1).*Q1c);
            CbAC(a3,a4,:) = zigzag(A21);
            A22=round(A12(x1,y1).*Q1c);
            CrAC(a4,a4,:) = zigzag(A22);
        end
    end
end
% ZRL 
YZRL = cell(M1,N1);
CbZRL = cell(M2,N2);
CrZRL = cell(M2,N2);
for a1 = 1:M1
    for a2 = 1:N1
        YZRL{a1,a2} = ZRLenco(YAC(a1,a2,:));
        if (a1<=M2)&&(a2<=N2)
            CbZRL{a1,a2} = ZRLenco(CbAC(a1,a2,:));
            CrZRL{a1,a2} = ZRLenco(CrAC(a1,a2,:));
        end
    end
end
% Process ZRL into (Run/group/index)
% set up group table
group(1,1:2)=[-1 1];
group(2,1:4)=[-3:-2 2:3];
group(3,1:8)=[-7:-4 4:7];
group(4,1:16)=[-15:-8 8:15];
group(5,1:32)=[-31:-16 16:31];
group(6,1:64)=[-63:-32 32:63];
group(7,1:128)=[-127:-64 64:127];
group(8,1:256)=[-255:-128 128:255];
group(9,1:512)=[-511:-256 256:511];
group(10,1:1024)=[-1023:-512 512:1023];
group(11,1:2048)=[-2047:-1024 1024:2047];
YZRL2 = cell(M1,N1);
CbZRL2 = cell(M2,N2);
CrZRL2 = cell(M2,N2);
for a1 = 1:M1
    for a2 = 1:N1
        YZRL2{a1,a2} = zeros(3,size(YZRL{a1,a2},2));
        YZRL2{a1,a2}(1,:) = YZRL{a1,a2}(1,:); % the row containing the zero run is the same
        for a3 = 1:size(YZRL{a1,a2},2)
            if YZRL{a1,a2}(2,a3) ~= 0
                [grp id] = find(group==YZRL{a1,a2}(2,a3));
                YZRL2{a1,a2}(2,a3) = grp;
                YZRL2{a1,a2}(3,a3) = id;
            end
        end
        if (a1<=M2)&&(a2<=N2)
            CbZRL2{a1,a2} = zeros(3,size(CbZRL{a1,a2},2));
            CbZRL2{a1,a2}(1,:) = CbZRL{a1,a2}(1,:);
            for a3 = 1:size(CbZRL{a1,a2},2)
                if CbZRL{a1,a2}(2,a3) ~= 0
                    [grp id] = find(group==CbZRL{a1,a2}(2,a3));
                    CbZRL2{a1,a2}(2,a3) = grp;
                    CbZRL2{a1,a2}(3,a3) = id;
                end
            end
            CrZRL2{a1,a2} = zeros(3,size(CrZRL{a1,a2},2));
            CrZRL2{a1,a2}(1,:) = CrZRL{a1,a2}(1,:);
            for a3 = 1:size(CrZRL{a1,a2},2)
                if CrZRL{a1,a2}(2,a3) ~= 0
                    [grp id] = find(group==CrZRL{a1,a2}(2,a3));
                    CrZRL2{a1,a2}(2,a3) = grp;
                    CrZRL2{a1,a2}(3,a3) = id;
                end
            end
        end
    end
end
% frequency table setup
% symbols:(run, group index)          : (0,0) (15,0) (0,1) (0,2) ... (0,10) (1,1)
% mapping each symbol to 1-D array    :   1      2     3     4   ...   12     13
%               mapping relation      : (a,b) -> (a*maxgroup+b+2)   *except for (a,b) = (0,0),(15,0)
%                freq table size      : 16*maxgroup + 2
[maxYgrp dummy] = find(group==max(max(max(abs(YAC)))));
[maxCbgrp dummy] = find(group==max(max(max(abs(CbAC)))));
[maxCrgrp dummy] = find(group==max(max(max(abs(CrAC)))));
Ycounts = ones(1,16*maxYgrp+2); 
Cbcounts = ones(1,16*maxCbgrp+2);
Crcounts = ones(1,16*maxCrgrp+2);
Ygrpidcounts = zeros(maxYgrp,2^maxYgrp);
Cbgrpidcounts = zeros(maxCbgrp,2^maxCbgrp);
Crgrpidcounts = zeros(maxCrgrp,2^maxCrgrp);
% initialize freq table
Ycounts(1) = 20;    Cbcounts(1) = 20;   Crcounts(1) = 20;
% run/grp table
for run = 0:15
    for grp = 1:maxYgrp
        id = run*maxYgrp + grp + 2;
        Ycounts(id) = 50/(2^(run+grp));
    end
    for grp = 1:maxCbgrp
        id = run*maxCbgrp + grp + 2;
        Cbcounts(id) = 5/(2^(run+grp));
    end
    for grp = 1:maxCrgrp
        id = run*maxCrgrp + grp + 2;
        Crcounts(id) = 5/(2^(run+grp));
    end
end
% grp index table
scale = 10;
dev = 1.2;
for i = 1:maxYgrp
    Ygrpidcounts(i,1:2^i) = ceil(scale*fspecial('gaussian',[1,2^i],dev)) + 1;
end
for i = 1:maxCbgrp
    Cbgrpidcounts(i,1:2^i) = ceil(scale*fspecial('gaussian',[1,2^i],dev)) + 1;
end
for i = 1:maxCrgrp
    Crgrpidcounts(i,1:2^i) = ceil(scale*fspecial('gaussian',[1,2^i],dev)) + 1;
end
% Parameters
G = 1.0006; % growth factor
code = zeros(1,M*N*8*3);
code_loc = 0;
data_high = 1;
data_low = 0;
ADD = ones(1,2); % [YADD CADD]
group_bits = 0;
% Main
for i = 1:M1
    for j = 1:N1
        num = 0;
        flag150 = 0;
        for k = 1:size(YZRL2{i,j},2)
            run = YZRL2{i,j}(1,k);
            grp = YZRL2{i,j}(2,k);
            grp_id = YZRL2{i,j}(3,k);
            if run==0 && grp==0
                id = 1;
            elseif run==15 && grp==0
                id = 2;
            else
                id = run*maxYgrp + grp + 2;
            end
            if flag150 == 0
                local_table = Ycounts;
                local_table(1) = local_table(1)*(num/16+1);
                s = cumsum(local_table/sum(local_table));
            else
                local_table = Ycounts;
                local_table(1) = 0; % if the previous ZRL pair is (15,0), the next pair is guaranteed not to be EOB (0,0)
                s = cumsum(local_table/sum(local_table));
                flag150 = 0;% restore the flag
            end
            num = num + run + 1;
            if (YZRL2{i,j}(1,k)==15) && (YZRL2{i,j}(2,k)==0)
               flag150 = 1;
            end
            % Range adjusting for run/grp: for example while encoding (1,1), we increase (2,1), (3,1) etc.
            % exception: 0/0  15/0
            if (id~=1) && (id~=2)
                for run2 = 0:15
                    id2 = run2*maxYgrp + grp + 2;
                    dist = abs(run2 - run) + 1;
                    if run2 ~= run
                        Ycounts(id2) = Ycounts(id2) + ADD(1)/(dist^6);
                    end
                end
            end
            % Range adjusting for run/grp: for example while encoding (1,1), we increase (1,2), (1,3) etc.
            % exception: 0/0  15/0
            if (id~=1) && (id~=2)
                for grp2 = 0:maxYgrp
                    id2 = run*maxYgrp + grp2 + 2;
                    dist = abs(grp2 - grp) + 1;
                    if grp2 ~= grp
                        Ycounts(id2) = Ycounts(id2) + ADD(1)/(dist^9);
                    end
                end
            end
            Ycounts(id) = Ycounts(id) + ADD(1);
            prevupper = data_high;
            prevlower = data_low;
            if id ~= 1
                data_low = data_low+(prevupper-prevlower)*s(id-1);
            end
            data_high = prevlower+(prevupper-prevlower)*s(id);
            % Rescaling 
            while (data_high<=0.5)||(data_low>=0.5)
                  if data_high<=0.5
                       code_loc = code_loc + 1;
                       code(code_loc) = 0;
                       data_high = 2 * data_high;
                       data_low =  2 * data_low;
                  else
                       code_loc = code_loc + 1;
                       code(code_loc) = 1;
                       data_high = 2 * data_high -1;
                       data_low =  2 * data_low -1;
                  end
                  if data_high == data_low
                      fprintf('ERROR: upper = lower occurs !\n');
                      return;
                  end
            end
            %encode Y group id if grp>0
            if grp>0
                s = cumsum(Ygrpidcounts(grp,:)/sum(Ygrpidcounts(grp,:)));
                for grp_id2 = 1:2^grp
                    dist = abs(grp_id2 - grp_id) + 1;
                    Ygrpidcounts(grp,grp_id2) = Ygrpidcounts(grp,grp_id2) + ADD(1)/(dist^3);
                end
                prevupper = data_high;
                prevlower = data_low;
                if grp_id ~= 1
                    data_low = data_low+(prevupper-prevlower)*s(grp_id-1);
                end
                data_high = prevlower+(prevupper-prevlower)*s(grp_id);
                while (data_high<=0.5)||(data_low>=0.5)
                    if data_high<=0.5
                        code_loc = code_loc + 1;
                        code(code_loc) = 0;
                        data_high = 2 * data_high;
                        data_low =  2 * data_low;
                    else
                        code_loc = code_loc + 1;
                        code(code_loc) = 1;
                        data_high = 2 * data_high -1;
                        data_low =  2 * data_low -1;
                    end
                    if data_high == data_low
                        fprintf('ERROR: upper = lower occurs !\n');
                        return;
                    end
                end
            end
            ADD(1) = ADD(1) * G;
            if ADD(1) > 1000
                Ycounts = ceil(Ycounts/2);
                ADD(1) = ADD(1)/2;
            end
        end
        if (i<=M2)&&(j<=N2)
            num = 0;
            flag150 = 0;
            for k = 1:size(CbZRL2{i,j},2)
                run = CbZRL2{i,j}(1,k);
                grp = CbZRL2{i,j}(2,k);
                grp_id = CbZRL2{i,j}(3,k);
                if run==0 && grp==0
                    id = 1;
                elseif run==15 && grp==0
                    id = 2;
                else
                    id = run*maxCbgrp + grp + 2;
                end
                if flag150 == 0
                    local_table = Cbcounts;
                    local_table(1) = local_table(1)*(num/16+1);
                    s = cumsum(local_table/sum(local_table));
                else
                    local_table = Cbcounts;
                    local_table(1) = 0; % if the previous ZRL pair is (15,0), the next pair is guaranteed not to be EOB (0,0)
                    s = cumsum(local_table/sum(local_table));
                    flag150 = 0;% restore the flag
                end
                num = num + run + 1;
                if (CbZRL2{i,j}(1,k)==15) && (CbZRL2{i,j}(2,k)==0)
                    flag150 = 1;
                end
                if (id~=1) && (id~=2)
                    for run2 = 0:15
                        id2 = run2*maxCbgrp + grp + 2;
                        dist = abs(run2 - run) + 1;
                        if run2 ~= run
                            Cbcounts(id2) = Cbcounts(id2) + ADD(2)/(dist^6);
                        end
                    end
                end
                if (id~=1) && (id~=2)
                    for grp2 = 0:maxCbgrp
                        id2 = run*maxCbgrp + grp2 + 2;
                        dist = abs(grp2 - grp) + 1;
                        if grp2 ~= grp
                            Cbcounts(id2) = Cbcounts(id2) + ADD(2)/(dist^9);
                        end
                    end
                end
                Cbcounts(id) = Cbcounts(id) + ADD(2);
                prevupper = data_high;
                prevlower = data_low;
                if id ~= 1
                    data_low = data_low+(prevupper-prevlower)*s(id-1);
                end
                data_high = prevlower+(prevupper-prevlower)*s(id);
                % Rescaling 
                while (data_high<=0.5)||(data_low>=0.5)
                    if data_high<=0.5
                        code_loc = code_loc + 1;
                        code(code_loc) = 0;
                        data_high = 2 * data_high;
                        data_low =  2 * data_low;
                    else
                        code_loc = code_loc + 1;
                        code(code_loc) = 1;
                        data_high = 2 * data_high -1;
                        data_low =  2 * data_low -1;
                    end
                    if data_high == data_low
                        fprintf('ERROR: upper = lower occurs !\n');
                        return;
                    end
                end
                %encode Cb group id if grp>0
                if grp>0
                    s = cumsum(Cbgrpidcounts(grp,:)/sum(Cbgrpidcounts(grp,:)));
                    for grp_id2 = 1:2^grp
                        dist = abs(grp_id2 - grp_id) + 1;
                        Cbgrpidcounts(grp,grp_id2) = Cbgrpidcounts(grp,grp_id2) + ADD(2)/(dist^3);
                    end
                    if grp <= maxCrgrp
                        for grp_id2 = 1:2^grp
                            dist = abs(grp_id2 - grp_id) + 1;
                            Crgrpidcounts(grp,grp_id2) = Crgrpidcounts(grp,grp_id2) + ADD(2)/(dist^5);
                        end
                    end
                    prevupper = data_high;
                    prevlower = data_low;
                    if grp_id ~= 1
                        data_low = data_low+(prevupper-prevlower)*s(grp_id-1);
                    end
                    data_high = prevlower+(prevupper-prevlower)*s(grp_id);
                    while (data_high<=0.5)||(data_low>=0.5)
                        if data_high<=0.5
                            code_loc = code_loc + 1;
                            code(code_loc) = 0;
                            data_high = 2 * data_high;
                            data_low =  2 * data_low;
                        else
                            code_loc = code_loc + 1;
                            code(code_loc) = 1;
                            data_high = 2 * data_high -1;
                            data_low =  2 * data_low -1;
                        end
                        if data_high == data_low
                            fprintf('ERROR: upper = lower occurs !\n');
                            return;
                        end
                    end
                end
                ADD(2) = ADD(2) * G;
                if ADD(2) > 100
                    Cbcounts = ceil(Cbcounts/2);
                    ADD(2) = ADD(2)/2;
                end
            end
            num = 0;
            flag150 = 0;
            for k = 1:size(CrZRL2{i,j},2)
                run = CrZRL2{i,j}(1,k);
                grp = CrZRL2{i,j}(2,k);
                grp_id = CrZRL2{i,j}(3,k);
                if run==0 && grp==0
                    id = 1;
                elseif run==15 && grp==0
                    id = 2;
                else
                    id = run*maxCrgrp + grp + 2;
                end
                if flag150 == 0
                    local_table = Crcounts;
                    local_table(1) = local_table(1)*(num/16+1);
                    s = cumsum(local_table/sum(local_table));
                else
                    local_table = Crcounts;
                    local_table(1) = 0; % if the previous ZRL pair is (15,0), the next pair is guaranteed not to be EOB (0,0)
                    s = cumsum(local_table/sum(local_table));
                    flag150 = 0;% restore the flag
                end
                num = num + run + 1;
                if (CrZRL2{i,j}(1,k)==15) && (CrZRL2{i,j}(2,k)==0)
                    flag150 = 1;
                end
                if (id~=1) && (id~=2)
                    for run2 = 0:15
                        id2 = run2*maxCrgrp + grp + 2;
                        dist = abs(run2 - run) + 1;
                        if run2 ~= run
                            Crcounts(id2) = Crcounts(id2) + ADD(2)/(dist^6);
                        end
                    end
                end
                if (id~=1) && (id~=2)
                    for grp2 = 0:maxCrgrp
                        id2 = run*maxCrgrp + grp2 + 2;
                        dist = abs(grp2 - grp) + 1;
                        if grp2 ~= grp
                            Crcounts(id2) = Crcounts(id2) + ADD(2)/(dist^9);
                        end
                    end
                end
                Crcounts(id) = Crcounts(id) + ADD(2);
                prevupper = data_high;
                prevlower = data_low;
                if id ~= 1
                    data_low = data_low+(prevupper-prevlower)*s(id-1);
                end
                data_high = prevlower+(prevupper-prevlower)*s(id);
                % Rescaling 
                while (data_high<=0.5)||(data_low>=0.5)
                    if data_high<=0.5
                        code_loc = code_loc + 1;
                        code(code_loc) = 0;
                        data_high = 2 * data_high;
                        data_low =  2 * data_low;
                    else
                        code_loc = code_loc + 1;
                        code(code_loc) = 1;
                        data_high = 2 * data_high -1;
                        data_low =  2 * data_low -1;
                    end
                    if data_high == data_low
                        fprintf('ERROR: upper = lower occurs !\n');
                        return;
                    end
                end
                %encode Cr group id if grp>0
                if grp>0
                    %group_bits = group_bits + grp;
                    s = cumsum(Crgrpidcounts(grp,:)/sum(Crgrpidcounts(grp,:)));
                    for grp_id2 = 1:2^grp
                        dist = abs(grp_id2 - grp_id) + 1;
                        Crgrpidcounts(grp,grp_id2) = Crgrpidcounts(grp,grp_id2) + ADD(2)/(dist^3);
                    end
                    if grp <= maxCbgrp
                        for grp_id2 = 1:2^grp
                            dist = abs(grp_id2 - grp_id) + 1;
                            Cbgrpidcounts(grp,grp_id2) = Cbgrpidcounts(grp,grp_id2) + ADD(2)/(dist^5);
                        end
                    end
                    prevupper = data_high;
                    prevlower = data_low;
                    if grp_id ~= 1
                        data_low = data_low+(prevupper-prevlower)*s(grp_id-1);
                    end
                    data_high = prevlower+(prevupper-prevlower)*s(grp_id);
                    while (data_high<=0.5)||(data_low>=0.5)
                        if data_high<=0.5
                            code_loc = code_loc + 1;
                            code(code_loc) = 0;
                            data_high = 2 * data_high;
                            data_low =  2 * data_low;
                        else
                            code_loc = code_loc + 1;
                            code(code_loc) = 1;
                            data_high = 2 * data_high -1;
                            data_low =  2 * data_low -1;
                        end
                        if data_high == data_low
                            fprintf('ERROR: upper = lower occurs !\n');
                            return;
                        end
                    end
                end
                ADD(2) = ADD(2) * G;
                if ADD(2) > 100
                    Crcounts = ceil(Crcounts/2);
                    ADD(2) = ADD(2)/2;
                end
            end
        end
    end
end
% binary encoding
code(code_loc+1:code_loc+2) = [0 1];
code_loc = code_loc + 2;
% truncate
code = code(1:code_loc);
info_bits = 4*3;% for maxY,Cb,Crgroup
fprintf('Proposed Method : code length = %g bits\n',length(code)+info_bits);
result(img) = length(code)+info_bits;
toc

% Decoder ---------------------------------------------------------------------------------------------------------------------------------
% Frequency table setup
Ycounts = ones(1,16*maxYgrp+2); 
Cbcounts = ones(1,16*maxCbgrp+2);
Crcounts = ones(1,16*maxCrgrp+2);
Ygrpidcounts = zeros(maxYgrp,2^maxYgrp);
Cbgrpidcounts = zeros(maxCbgrp,2^maxCbgrp);
Crgrpidcounts = zeros(maxCrgrp,2^maxCrgrp);
% initialize freq table
Ycounts(1) = 20;    Cbcounts(1) = 20;   Crcounts(1) = 20;
% run/grp table
for run = 0:15
    for grp = 1:maxYgrp
        id = run*maxYgrp + grp + 2;
        Ycounts(id) = 50/(2^(run+grp));
    end
    for grp = 1:maxCbgrp
        id = run*maxCbgrp + grp + 2;
        Cbcounts(id) = 5/(2^(run+grp));
    end
    for grp = 1:maxCrgrp
        id = run*maxCrgrp + grp + 2;
        Crcounts(id) = 5/(2^(run+grp));
    end
end
% grp index table
scale = 10;
dev = 1.2;
for i = 1:maxYgrp
    Ygrpidcounts(i,1:2^i) = ceil(scale*fspecial('gaussian',[1,2^i],dev)) + 1;
end
for i = 1:maxCbgrp
    Cbgrpidcounts(i,1:2^i) = ceil(scale*fspecial('gaussian',[1,2^i],dev)) + 1;
end
for i = 1:maxCrgrp
    Crgrpidcounts(i,1:2^i) = ceil(scale*fspecial('gaussian',[1,2^i],dev)) + 1;
end
% Parameters
G = 1.0006; % growth factor
ADD = ones(1,2); % [YADD CADD]
data_high = 1; %the temp of upper bound
data_low  = 0; %the temp of lower bound
bit_high  = 1; %the second temp of upper bound
bit_low   = 0; %the second temp of lower bound
code_loc=0; % record the location of code
bit_scale=0.5; % scaling number
YZRL2r = cell(M1,N1);
CbZRL2r = cell(M2,N2);
CrZRL2r = cell(M2,N2);
% Main
for i = 1:M1
    for j = 1:N1
        num = 0;
        flag150 = 0;
        for k = 1:size(YZRL2{i,j},2)
            flag=1;
            while(flag)% decoding run/grp
                if flag150 == 0
                    local_table = Ycounts;
                    local_table(1) = local_table(1)*(num/16+1);
                    s = cumsum(local_table/sum(local_table));
                else
                    local_table = Ycounts;
                    local_table(1) = 0; % if the previous ZRL pair is (15,0), the next pair is guaranteed not to be EOB (0,0)
                    s = cumsum(local_table/sum(local_table));
                end
                prevupper = data_high;
                prevlower = data_low;
                for ii=1:size(Ycounts,2)
                    templower=data_low;
                    if ii~=1
                        templower=templower+(prevupper-prevlower)*s(ii-1);
                    end
                    tempupper=prevlower+(prevupper-prevlower)*s(ii);
                    if templower<=bit_low && bit_high<=tempupper
                        id=ii;
                        flag=0;
                        data_high=tempupper;
                        data_low=templower;
                    end
                end
                if(flag == 0)
                    break;
                end
                code_loc=code_loc+1;
                if code_loc>length(code)
                    code(code_loc)=1;
                end
                if code(code_loc)==1
                    bit_low  = bit_low+bit_scale;
                else
                    bit_high = bit_high-bit_scale;
                end
                bit_scale=bit_scale*0.5;
            end
            % post processing after decoding run/grp
            if flag150 == 1
                flag150 = 0;
            end
            if id == 1
                run = 0;
                grp = 0;
            elseif id == 2
                run = 15;
                grp = 0;
            else
                % id = run*maxYgrp + grp + 2  (grp: 1~maxYgrp)
                for run2 = 0:15
                    for grp2 = 1:maxYgrp
                        id2 = run2*maxYgrp + grp2 + 2;
                        if id2 == id
                            run = run2;
                            grp = grp2;
                        end
                    end
                end
            end
            YZRL2r{i,j}(1,k) = run;
            YZRL2r{i,j}(2,k) = grp;
            num = num + run + 1;
            if (YZRL2r{i,j}(1,k)==15) && (YZRL2r{i,j}(2,k)==0)
               flag150 = 1;
            end
            % Range adjusting
            if (id~=1) && (id~=2)
                for run2 = 0:15
                    id2 = run2*maxYgrp + grp + 2;
                    dist = abs(run2 - run) + 1;
                    if run2 ~= run
                        Ycounts(id2) = Ycounts(id2) + ADD(1)/(dist^6);
                    end
                end
            end
            if (id~=1) && (id~=2)
                for grp2 = 0:maxYgrp
                    id2 = run*maxYgrp + grp2 + 2;
                    dist = abs(grp2 - grp) + 1;
                    if grp2 ~= grp
                        Ycounts(id2) = Ycounts(id2) + ADD(1)/(dist^9);
                    end
                end
            end
            Ycounts(id) = Ycounts(id) + ADD(1);
            while (data_high<=0.5)||(data_low>=0.5)
                if data_high<=0.5
                    data_high=2*data_high;
                    data_low=2*data_low;
                    bit_high=2*bit_high;
                    bit_low=2*bit_low;
                else
                    data_high = 2*data_high -1;
                    data_low  = 2*data_low -1;
                    bit_high  = 2*bit_high-1;
                    bit_low   = 2*bit_low-1;
                end
                bit_scale = 2*bit_scale;   %bit backward
                if data_high == data_low
                    fprintf('ERROR: upper = lower occurs !\n');
                    return;
                end
            end
            % decode Y group id if grp>0
            grp_id = 0;% for cases like 0/0 15/0
            if grp>0
                flag = 1;
                while(flag)
                    s = cumsum(Ygrpidcounts(grp,:)/sum(Ygrpidcounts(grp,:)));
                    prevupper = data_high;
                    prevlower = data_low;
                    for ii=1:2^grp
                        templower=data_low;
                        if ii~=1
                            templower=templower+(prevupper-prevlower)*s(ii-1);
                        end
                        tempupper=prevlower+(prevupper-prevlower)*s(ii);
                        if templower<=bit_low && bit_high<=tempupper
                            grp_id=ii;
                            flag=0;
                            data_high=tempupper;
                            data_low=templower;
                        end
                    end
                    if(flag == 0)
                        break;
                    end
                    code_loc=code_loc+1;
                    if code_loc>length(code)
                        code(code_loc)=1;
                    end
                    if code(code_loc)==1
                        bit_low  = bit_low+bit_scale;
                    else
                        bit_high = bit_high-bit_scale;
                    end
                    bit_scale=bit_scale*0.5;
                end
                % post processing after decoding grp_id
                for grp_id2 = 1:2^grp
                    dist = abs(grp_id2 - grp_id) + 1;
                    Ygrpidcounts(grp,grp_id2) = Ygrpidcounts(grp,grp_id2) + ADD(1)/(dist^3);
                end
                while (data_high<=0.5)||(data_low>=0.5)
                    if data_high<=0.5
                        data_high=2*data_high;
                        data_low=2*data_low;
                        bit_high=2*bit_high;
                        bit_low=2*bit_low;
                    else
                        data_high = 2*data_high -1;
                        data_low  = 2*data_low -1;
                        bit_high  = 2*bit_high-1;
                        bit_low   = 2*bit_low-1;
                    end
                    bit_scale = 2*bit_scale;   %bit backward
                    if data_high == data_low
                        fprintf('ERROR: upper = lower occurs !\n');
                        return;
                    end
                end
            end
            YZRL2r{i,j}(3,k) = grp_id;
            ADD(1) = ADD(1) * G;
            if ADD(1) > 1000
                Ycounts = ceil(Ycounts/2);
                ADD(1) = ADD(1)/2;
            end
        end
        if (i<=M2)&&(j<=N2)
            num = 0;
            flag150 = 0;
            for k = 1:size(CbZRL2{i,j},2)
                flag=1;
                while(flag)% decoding run/grp
                    if flag150 == 0
                        local_table = Cbcounts;
                        local_table(1) = local_table(1)*(num/16+1);
                        s = cumsum(local_table/sum(local_table));
                    else
                        local_table = Cbcounts;
                        local_table(1) = 0; % if the previous ZRL pair is (15,0), the next pair is guaranteed not to be EOB (0,0)
                        s = cumsum(local_table/sum(local_table));
                    end
                    prevupper = data_high;
                    prevlower = data_low;
                    for ii=1:size(Cbcounts,2)
                        templower=data_low;
                        if ii~=1
                            templower=templower+(prevupper-prevlower)*s(ii-1);
                        end
                        tempupper=prevlower+(prevupper-prevlower)*s(ii);
                        if templower<=bit_low && bit_high<=tempupper
                            id=ii;
                            flag=0;
                            data_high=tempupper;
                            data_low=templower;
                        end
                    end
                    if(flag == 0)
                        break;
                    end
                    code_loc=code_loc+1;
                    if code_loc>length(code)
                        code(code_loc)=1;
                    end
                    if code(code_loc)==1
                        bit_low  = bit_low+bit_scale;
                    else
                        bit_high = bit_high-bit_scale;
                    end
                    bit_scale=bit_scale*0.5;
                end
                % post processing after decoding run/grp
                if flag150 == 1
                    flag150 = 0;
                end
                if id == 1
                    run = 0;
                    grp = 0;
                elseif id == 2
                    run = 15;
                    grp = 0;
                else
                    for run2 = 0:15
                        for grp2 = 1:maxCbgrp
                            id2 = run2*maxCbgrp + grp2 + 2;
                            if id2 == id
                                run = run2;
                                grp = grp2;
                            end
                        end
                    end
                end
                CbZRL2r{i,j}(1,k) = run;
                CbZRL2r{i,j}(2,k) = grp;
                num = num + run + 1;
                if (CbZRL2r{i,j}(1,k)==15) && (CbZRL2r{i,j}(2,k)==0)
                    flag150 = 1;
                end
                % Range adjusting
                if (id~=1) && (id~=2)
                    for run2 = 0:15
                        id2 = run2*maxCbgrp + grp + 2;
                        dist = abs(run2 - run) + 1;
                        if run2 ~= run
                            Cbcounts(id2) = Cbcounts(id2) + ADD(2)/(dist^6);
                        end
                    end
                end
                if (id~=1) && (id~=2)
                    for grp2 = 0:maxCbgrp
                        id2 = run*maxCbgrp + grp2 + 2;
                        dist = abs(grp2 - grp) + 1;
                        if grp2 ~= grp
                            Cbcounts(id2) = Cbcounts(id2) + ADD(2)/(dist^9);
                        end
                    end
                end
                Cbcounts(id) = Cbcounts(id) + ADD(2);
                while (data_high<=0.5)||(data_low>=0.5)
                    if data_high<=0.5
                        data_high=2*data_high;
                        data_low=2*data_low;
                        bit_high=2*bit_high;
                        bit_low=2*bit_low;
                    else
                        data_high = 2*data_high -1;
                        data_low  = 2*data_low -1;
                        bit_high  = 2*bit_high-1;
                        bit_low   = 2*bit_low-1;
                    end
                    bit_scale = 2*bit_scale;   %bit backward
                    if data_high == data_low
                        fprintf('ERROR: upper = lower occurs !\n');
                        return;
                    end
                end
                % decode Cb group id if grp>0
                grp_id = 0;% for cases like 0/0 15/0
                if grp>0
                    flag = 1;
                    while(flag)
                        s = cumsum(Cbgrpidcounts(grp,:)/sum(Cbgrpidcounts(grp,:)));
                        prevupper = data_high;
                        prevlower = data_low;
                        for ii=1:2^grp
                            templower=data_low;
                            if ii~=1
                                templower=templower+(prevupper-prevlower)*s(ii-1);
                            end
                            tempupper=prevlower+(prevupper-prevlower)*s(ii);
                            if templower<=bit_low && bit_high<=tempupper
                                grp_id=ii;
                                flag=0;
                                data_high=tempupper;
                                data_low=templower;
                            end
                        end
                        if(flag == 0)
                            break;
                        end
                        code_loc=code_loc+1;
                        if code_loc>length(code)
                            code(code_loc)=1;
                        end
                        if code(code_loc)==1
                            bit_low  = bit_low+bit_scale;
                        else
                            bit_high = bit_high-bit_scale;
                        end
                        bit_scale=bit_scale*0.5;
                    end
                    % post processing after decoding grp_id
                    for grp_id2 = 1:2^grp
                        dist = abs(grp_id2 - grp_id) + 1;
                        Cbgrpidcounts(grp,grp_id2) = Cbgrpidcounts(grp,grp_id2) + ADD(2)/(dist^3);
                    end
                    if grp <= maxCrgrp
                        for grp_id2 = 1:2^grp
                            dist = abs(grp_id2 - grp_id) + 1;
                            Crgrpidcounts(grp,grp_id2) = Crgrpidcounts(grp,grp_id2) + ADD(2)/(dist^5);
                        end
                    end
                    while (data_high<=0.5)||(data_low>=0.5)
                        if data_high<=0.5
                            data_high=2*data_high;
                            data_low=2*data_low;
                            bit_high=2*bit_high;
                            bit_low=2*bit_low;
                        else
                            data_high = 2*data_high -1;
                            data_low  = 2*data_low -1;
                            bit_high  = 2*bit_high-1;
                            bit_low   = 2*bit_low-1;
                        end
                        bit_scale = 2*bit_scale;   %bit backward
                        if data_high == data_low
                            fprintf('ERROR: upper = lower occurs !\n');
                            return;
                        end
                    end
                end
                CbZRL2r{i,j}(3,k) = grp_id;
                ADD(2) = ADD(2) * G;
                if ADD(2) > 100
                    Cbcounts = ceil(Cbcounts/2);
                    ADD(2) = ADD(2)/2;
                end
            end
            num = 0;
            flag150 = 0;
            for k = 1:size(CrZRL2{i,j},2)
                flag=1;
                while(flag)% decoding run/grp
                    if flag150 == 0
                        local_table = Crcounts;
                        local_table(1) = local_table(1)*(num/16+1);
                        s = cumsum(local_table/sum(local_table));
                    else
                        local_table = Crcounts;
                        local_table(1) = 0; % if the previous ZRL pair is (15,0), the next pair is guaranteed not to be EOB (0,0)
                        s = cumsum(local_table/sum(local_table));
                    end
                    prevupper = data_high;
                    prevlower = data_low;
                    for ii=1:size(Crcounts,2)
                        templower=data_low;
                        if ii~=1
                            templower=templower+(prevupper-prevlower)*s(ii-1);
                        end
                        tempupper=prevlower+(prevupper-prevlower)*s(ii);
                        if templower<=bit_low && bit_high<=tempupper
                            id=ii;
                            flag=0;
                            data_high=tempupper;
                            data_low=templower;
                        end
                    end
                    if(flag == 0)
                        break;
                    end
                    code_loc=code_loc+1;
                    if code_loc>length(code)
                        code(code_loc)=1;
                    end
                    if code(code_loc)==1
                        bit_low  = bit_low+bit_scale;
                    else
                        bit_high = bit_high-bit_scale;
                    end
                    bit_scale=bit_scale*0.5;
                end
                % post processing after decoding run/grp
                if flag150 == 1
                    flag150 = 0;
                end
                if id == 1
                    run = 0;
                    grp = 0;
                elseif id == 2
                    run = 15;
                    grp = 0;
                else
                    for run2 = 0:15
                        for grp2 = 1:maxCrgrp
                            id2 = run2*maxCrgrp + grp2 + 2;
                            if id2 == id
                                run = run2;
                                grp = grp2;
                            end
                        end
                    end
                end
                CrZRL2r{i,j}(1,k) = run;
                CrZRL2r{i,j}(2,k) = grp;
                num = num + run + 1;
                if (CrZRL2r{i,j}(1,k)==15) && (CrZRL2r{i,j}(2,k)==0)
                    flag150 = 1;
                end
                % Range adjusting
                if (id~=1) && (id~=2)
                    for run2 = 0:15
                        id2 = run2*maxCrgrp + grp + 2;
                        dist = abs(run2 - run) + 1;
                        if run2 ~= run
                            Crcounts(id2) = Crcounts(id2) + ADD(2)/(dist^6);
                        end
                    end
                end
                if (id~=1) && (id~=2)
                    for grp2 = 0:maxCrgrp
                        id2 = run*maxCrgrp + grp2 + 2;
                        dist = abs(grp2 - grp) + 1;
                        if grp2 ~= grp
                            Crcounts(id2) = Crcounts(id2) + ADD(2)/(dist^9);
                        end
                    end
                end
                Crcounts(id) = Crcounts(id) + ADD(2);
                while (data_high<=0.5)||(data_low>=0.5)
                    if data_high<=0.5
                        data_high=2*data_high;
                        data_low=2*data_low;
                        bit_high=2*bit_high;
                        bit_low=2*bit_low;
                    else
                        data_high = 2*data_high -1;
                        data_low  = 2*data_low -1;
                        bit_high  = 2*bit_high-1;
                        bit_low   = 2*bit_low-1;
                    end
                    bit_scale = 2*bit_scale;   %bit backward
                    if data_high == data_low
                        fprintf('ERROR: upper = lower occurs !\n');
                        return;
                    end
                end
                % decode Cr group id if grp>0
                grp_id = 0;% for cases like 0/0 15/0
                if grp>0
                    flag = 1;
                    while(flag)
                        s = cumsum(Crgrpidcounts(grp,:)/sum(Crgrpidcounts(grp,:)));
                        prevupper = data_high;
                        prevlower = data_low;
                        for ii=1:2^grp
                            templower=data_low;
                            if ii~=1
                                templower=templower+(prevupper-prevlower)*s(ii-1);
                            end
                            tempupper=prevlower+(prevupper-prevlower)*s(ii);
                            if templower<=bit_low && bit_high<=tempupper
                                grp_id=ii;
                                flag=0;
                                data_high=tempupper;
                                data_low=templower;
                            end
                        end
                        if(flag == 0)
                            break;
                        end
                        code_loc=code_loc+1;
                        if code_loc>length(code)
                            code(code_loc)=1;
                        end
                        if code(code_loc)==1
                            bit_low  = bit_low+bit_scale;
                        else
                            bit_high = bit_high-bit_scale;
                        end
                        bit_scale=bit_scale*0.5;
                    end
                    % post processing after decoding grp_id
                    for grp_id2 = 1:2^grp
                        dist = abs(grp_id2 - grp_id) + 1;
                        Crgrpidcounts(grp,grp_id2) = Crgrpidcounts(grp,grp_id2) + ADD(2)/(dist^3);
                    end
                    if grp <= maxCbgrp
                        for grp_id2 = 1:2^grp
                            dist = abs(grp_id2 - grp_id) + 1;
                            Cbgrpidcounts(grp,grp_id2) = Cbgrpidcounts(grp,grp_id2) + ADD(2)/(dist^5);
                        end
                    end
                    while (data_high<=0.5)||(data_low>=0.5)
                        if data_high<=0.5
                            data_high=2*data_high;
                            data_low=2*data_low;
                            bit_high=2*bit_high;
                            bit_low=2*bit_low;
                        else
                            data_high = 2*data_high -1;
                            data_low  = 2*data_low -1;
                            bit_high  = 2*bit_high-1;
                            bit_low   = 2*bit_low-1;
                        end
                        bit_scale = 2*bit_scale;   %bit backward
                        if data_high == data_low
                            fprintf('ERROR: upper = lower occurs !\n');
                            return;
                        end
                    end
                end
                CrZRL2r{i,j}(3,k) = grp_id;
                ADD(2) = ADD(2) * G;
                if ADD(2) > 100
                    Crcounts = ceil(Crcounts/2);
                    ADD(2) = ADD(2)/2;
                end
            end
        end
    end
end
if ~isequal(YZRL2,YZRL2r) || ~isequal(CbZRL2,CbZRL2r) || ~isequal(CrZRL2,CrZRL2r)
    fprintf('image: %g invalid\n',img);
    break;
end

end
result