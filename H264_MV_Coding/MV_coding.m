% MV_Coding
% [Input]: the "Motion Vectors Differencs (MVDs)" with quarter-accuracy from 
% 13 testing videos (the first 100 frames), the motion vector is obtained 
% by the H.264/AVC Rate Distortion Optimization rule [31], with QP = 28  
clc; clear;
% Main
datalist_QCIF = {'Akiyo_QCIF','Carphone_QCIF','City_QCIF','Coastguard_QCIF','Grandma_QCIF','Hall_QCIF','Ice_QCIF','Suzie_QCIF'};
datalist_CIF = {'Foreman_CIF','MissAmerica_CIF','MobileCalendar_CIF','Stefan_CIF','Waterfall_CIF'};
addpath(genpath('.'));
Final_Result = cell(15,8);
Final_Result(1,:) = {[], 'SAC', [], 'AAC', [], 'Pro', [], 'Entropy'}; 
Final_Result(2,:) = {[], 'BITS PER MVD','SUCCESS in ENC/DEC?', 'BITS PER MVD','SUCCESS in ENC/DEC?', 'BITS PER MVD','SUCCESS in ENC/DEC?', 'E(X)+E(Y)'}; 

total_bit_sac = 0; % total bit per MV (SAC)
total_bit_aac = 0; % total bit per MV (AAC)
total_bit_pro = 0; % total bit per MV (Proposed)
total_E = 0; % total entropy

fprintf('Program running...\n');
for num = 1:13

    vx = []; vy = [];
    XX = []; YY = [];
    YY1 =[]; XX1=[];
    
    bitstream_SAC = '';
    bitstream_AAC = '';
    bitstream_Pro = '';
    
    if num <= 8 % QCIF
        dat_path = datalist_QCIF{num};
        Final_Result(num+2,1) = {dat_path};
    else % num>=8  CIF
        dat_path = datalist_CIF{num-8};
        Final_Result(num+2,1) = {dat_path};
    end
        filename1 = fullfile(dat_path, 'MVx.mat');
        filename2 = fullfile(dat_path, 'MVy.mat');
        load(filename1,'-mat','MVx'); % 1x149 Cell with double,quarter vector value
        load(filename2,'-mat','MVy'); % 1x149 Cell with double,quarter vector value
        
% (Pre) ===== SYMBOL PRE-PROCESSING =====
%  Input: MVDs with quarter-accuracy { Z + [ -0.25 ~ +0.75 ] }
% Output: Positive Integers {N}, the symbols for Proposed algorithm
    for i = 1:length(MVx)
        vx = [vx MVx{i} * 4]; % -0.25~-1 0.25~1
        vy = [vy MVy{i} * 4];
    end  
    parfor (j = 1:length(vx),2) % Parallel Computing by two Cores
        if vx(j) > 0
            XX(j) = 2 * vx(j);
        elseif vx(j) == 0
            XX(j) = 1;
        else % MVx(i,j) < 0
            XX(j) = (-2) * vx(j) +1;
        end
        if vy(j) > 0
            YY(j) = 2 * vy(j);
        elseif vx(j) == 0
            YY(j) = 1;
        else % MVx(i,j) < 0
            YY(j) = (-2) * vy(j) +1;
        end
    end
%% (1) =================== SAC ===================
    PX1 = histc(XX, 1:max(XX));
    PY1 = histc(YY, 1:max(YY));

    [bitstream_SAC,x_end] = arith(XX, PX1, 'x', 'SAC');
    bitstream_SAC = [bitstream_SAC arith(YY, PY1, 'y', 'SAC')];
    
    % Calculate compressed bitstream length (along with side info)
    bitlength_SAC = length(bitstream_SAC) + side_info_est_bits(PX1) + side_info_est_bits(PY1);
    total_bit_sac = total_bit_sac + bitlength_SAC/length(XX);
    
    XX_SAC = iarith(bitstream_SAC(1:x_end), PX1, length(XX), 'x', 'SAC');
    YY_SAC = iarith(bitstream_SAC(x_end+1:end), PY1, length(YY), 'y', 'SAC');
    
    Final_Result(num+2,2) = {bitlength_SAC/length(XX)}; % 'BITS PER MV'
    
    % Check Validity
    if (all(XX_SAC == XX)) && (all(YY_SAC == YY))
        Final_Result(num+2,3) = {'Success'};
    else
        Final_Result(num+2,3) = {'Mismatch'};
    end
    
%% (2) =================== AAC ===================
    PX2 = ones(1,max(XX));
    PY2 = ones(1,max(YY));
    
    [bitstream_AAC,x_end] = arith(XX, PX2, 'x', 'AAC');
    bitstream_AAC = [bitstream_AAC arith(YY,PY2, 'y', 'AAC')];
    total_bit_aac = total_bit_aac + length(bitstream_AAC)/length(XX);
    
    XX_AAC = iarith(bitstream_AAC(1:x_end), PX2, length(XX), 'x', 'AAC');
    YY_AAC = iarith(bitstream_AAC(x_end+1:end), PY2, length(YY), 'y', 'AAC');
    
    Final_Result(num+2,4) = {length(bitstream_AAC)/length(XX)}; % 'BITS PER MV'
    
    % Check Validity
    if (all(XX_AAC == XX)) && (all(YY_AAC == YY))
        Final_Result(num+2,5) = {'Success'};
    else
        Final_Result(num+2,5) = {'Mismatch'};
    end
    
%% (3) ===== Proposed MVD encoding =====
%  Input: MVDs information in Natural numbers, initial frequency table
% Output: bitstream constructed by {0,1}
    % ***** Section II. A. Initial Frequency Table *****
    PX = ini_freq(max(XX));
    PY = ini_freq(max(YY));
    
    [bitstream_Pro,x_end] = arith(XX, PX, 'x', 'Pro');
    bitstream_Pro = [bitstream_Pro arith(YY,PY, 'y', 'Pro')];
    
    Final_Result(num+2,6) = {length(bitstream_Pro)/length(XX)}; % 'BITS PER MV'
    total_bit_pro = total_bit_pro + length(bitstream_Pro)/length(XX);
% =================== Decoder ===================
% (2) ===== Proposed MVD decoding =====
%  Input: MVDs information in Natural numbers, initial frequency table
% Output: Positive Integers {N}, the symbols for Proposed algorithm
    PX = ini_freq(max(XX));
    PY = ini_freq(max(YY));
    
    XX_prop = iarith(bitstream_Pro(1:x_end),PX,length(XX),'x','Pro');
    YY_prop = iarith(bitstream_Pro(x_end+1:end),PY,length(YY),'y','Pro');
    
    % Check Validity
    if (all(XX_prop == XX)) && (all(YY_prop == YY))
        Final_Result(num+2,7) = {'Success'};
    else
        Final_Result(num+2,7) = {'Mismatch'};
    end
%% (4) ===== Entropy =====
    EX = calEntropy(XX,max(XX));
    EY = calEntropy(YY,max(YY));
    Final_Result(num+2,8) = {EX+EY};
    total_E = total_E + EX + EY;
end
Final_Result(num+3,1) = {'AVERAGE'};
Final_Result(num+3,2) = {total_bit_sac/num};
Final_Result(num+3,4) = {total_bit_aac/num};
Final_Result(num+3,6) = {total_bit_pro/num};
Final_Result(num+3,8) = {total_E/num};

fprintf('Done \n');

