clc;clear;
close all;
sign4online = 0;
method = {'DSB'                1;
          'DSB_Mask'           0;
          'MVDR'               0;
          'MVDR_ESB'           0;
          'MVDR_AESB'          0;
          'MVDR_PCA'           0;
          'MVDR_Search'        0;
          'MVDR_AESB_Search'   0;
          'LCMV'               0;
          'LCMV_ESB'           0;
          'LCMV_AESB'          0;
          'LCMV_Search'        0;
          'ICA_funda'          0;
          'ICA_initial'        0;
          'ICA_Sawada'         0;
          'IVA'                0;
          'IVA_woDR'           0;
          'AuxIVA'             0;
          'OverIVA'            0;
          'OverILRMA'          0;
          'ILRMA_woDR'         0;
          'ILRMA'              0;
          'ILRMA_PF'           0;
          'FastMNMF1'          0;
          'FastMNMF2'          0;
          'FastFCA_AS'         0; %%%% remaining to be written
          'maxSNR'             0;
          %%%%%%% Dereverberation
          'WPE'                0;
          %%%%%%% Compound method
          'cGMM_maxSNR'        0;
          'ALL'                0
          };                    %%%%%%% This version only support single method once, dont select multiple methods and 'ALL'
speakerNum = 10;
pathFe = '../simuData/sound sources/simulation/female';
dirFe = dir(pathFe);
pathMa = '../simuData/sound sources/simulation/male';
dirMa = dir(pathMa);
for T60 = 0:0.3:0.3%0.9
angle = 0;
subpath = strcat('T60_',num2str(T60),'/Angle_',num2str(angle));
sOri{speakerNum} = [];
sZero{speakerNum} = [];
headPathOri = '../simuData/sound sources/origin';
for i = 1:speakerNum
    if(i<=round(speakerNum/2))
        path = fullfile(pathFe,dirFe(i+2).name,subpath);
        pathOri = fullfile(headPathOri,'female',strcat(dirFe(i+2).name,'.wav'));
    else
        path = fullfile(pathMa,dirMa(i-speakerNum/2+2).name,subpath);
        pathOri = fullfile(headPathOri,'male',strcat(dirMa(i-speakerNum/2+2).name,'.wav'));
    end
    for j = 0:6
        signal(:,j+1) = audioread(strcat(path,'/Mic',num2str(j),'.wav'));
    end
    sZero{i} = signal;
    clear signal;
    [sOri{i},fs] = audioread(pathOri);       
end
MetricsAngle = {};
AngleNum = 1;
for angle = 0:5:180
    Metrics = [];
    sAngle{speakerNum} = [];
    subpath = strcat('T60_',num2str(T60),'/Angle_',num2str(angle));
    for i = 1:speakerNum
        if(i<=round(speakerNum/2))
            path = fullfile(pathFe,dirFe(i+2).name,subpath);
        else
            path = fullfile(pathMa,dirMa(i-speakerNum/2+2).name,subpath);
        end
        for j = 0:6
            signal(:,j+1) = audioread(strcat(path,'/Mic',num2str(j),'.wav'));
        end
        sAngle{i} = signal;
        clear signal;
    end
    %%%%%%%%%%%%%%%%%%%%%% Initialization
    SetupStruc = ISM_setup;
    SetupStruc.Angle = [0 angle];
    SetupStruc.fs = fs;
    for i = 1:speakerNum-1
        for j = i+1:speakerNum
            s1 = sZero{i};
            s2 = sAngle{j};
            sO1 = sOri{i};
            sO2 = sOri{j};
            Len = min([size(s1,1),size(s2,1),size(sO1,1),size(sO2,1)]);
            s1 = s1(1:Len,:);s2 = s2(1:Len,:);
            Eng1 = sqrt(sum(sum(s1.^2))/Len/SetupStruc.Channel_Num);
            Eng2 = sqrt(sum(sum(s2.^2))/Len/SetupStruc.Channel_Num);
            s1 = s1/Eng1; s2 = s2/Eng2;
            InSNR = 10*log10((sum(sum(s1.^2))/Len/SetupStruc.Channel_Num)/...
                (sum(sum(s2.^2))/Len/SetupStruc.Channel_Num));            
            sMix = s1+s2;
            S = [sO1(1:Len,:) sO2(1:Len,:)];
            if(sign4online==1)
                [Re,SetupStruc] = OnProcess(method,sMix,SetupStruc);
                eval(strcat('SE = Re.', Re.method{1}, '.OnS;'));
            else
                [Re,SetupStruc] = Process(method,sMix,SetupStruc);
                eval(strcat('SE = Re.', Re.method{1}, '.S;'));
            end
%             Unmix_s(:,:,1) = s1;
%             Unmix_s(:,:,2) = s2;
%             [Me] = Cal_metrics(Re,Unmix_s,[s1(:,1) s2(:,1)],SetupStruc);                    
            [SDR,SIR,SAR,perm]=bss_eval_sources(SE,S);
            Metrics = [Metrics;[SDR,SIR,SAR]];
        end
    end
    savePath = strcat('./Results/',Re.method{1},'/T60_',num2str(T60));
    if exist(savePath,'dir')==0
        mkdir(savePath);
    end
    eval(strcat('save(''./Results/',Re.method{1},'/T60_',num2str(T60),'/Metrics',num2str(angle),'.txt'',','''Metrics'',','''-ascii'');'));
    MetricsAngle{AngleNum} = Metrics;
    AngleNum = AngleNum+1;
    clear sAngle;
end
end









































