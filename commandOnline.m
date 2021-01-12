clc;clear;
close all;
Angle = [135,180];                               % Angle, Options:90,135,180,270 plus 0,45,225,315
[s,sOri,Unmix_s,SetupStruc] = readData(Angle,ISM_setup);    %'s' is the muli-channel mixture signal,size:Length*7,
                                                %'sOri' is the each single singal of the 1st channel corresponding to the angle,size:Length*N
method = {'DSB'                0;
          'DSB_Mask'           0;
          'MVDR'               0;
          'LCMV'               0;
          'maxSNR'             1;
          'ALL'                0
          };
SetupStruc.unS = Unmix_s;
[Re,SetupStruc] = OnProcess(method,s,SetupStruc);
[Me] = Cal_metrics(Re,Unmix_s,sOri,SetupStruc);