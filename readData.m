function [s,Ori_s,Unmix_s,SetupStruc] = readData(Angle,SetupStruc,channel)
Num = length(Angle);
SetupStruc.Angle = Angle;
T60 = SetupStruc.T60;
Sign_compared = SetupStruc.Sign_compared;
Channel_Num = SetupStruc.Channel_Num;
Len = zeros(Num,1);
if(~exist('channel','var'))
    channel = 1;
end
if(SetupStruc.Sign_Mid==1)
    channel_line = 0:Channel_Num-1;
else
    channel_line = 1:Channel_Num;
end
for n = 1:Num
    dataPath = 'simu_data/T60_';
    for i = channel_line%0:Channel_Num-1
        wavName =strcat(dataPath,num2str(T60),'/',num2str(Angle(n)),'/Mic',num2str(i),'.wav');
        if(SetupStruc.Sign_Mid==1)
            eval(strcat('[s_temp', num2str(n), '(:,', num2str(i), '+1),fs] = audioread(''', wavName, ''');'));
        else
            eval(strcat('[s_temp', num2str(n), '(:,', num2str(i), '),fs] = audioread(''', wavName, ''');'));
        end
    end
    eval(strcat('Len(n) = size(s_temp', num2str(n),',1);'));
end
SetupStruc.fs = fs;
len = min(Len);
s = s_temp1(1:len,:);
Unmix_s(:,:,1) = s;
if(Num>1)
    Ori_s = zeros(len,Num);   
    Ori_s(:,1) = s_temp1(1:len,channel);
else
    Ori_s = zeros(len,2);
    Ori_s_temp = audioread(strcat(dataPath,'0/',num2str(Angle(1)),'/Mic0.wav'));
    Ori_s(:,1) = Ori_s_temp(1:len);
    Ori_s(:,2) = s_temp1(1:len,channel);
end
for n = 2:Num
    eval(strcat('s_temp = s_temp', num2str(n), ';'));
    s = s+s_temp(1:len,:);
    Unmix_s(:,:,n) = s_temp(1:len,:);
    Ori_s(:,n) = s_temp(1:len,channel);
end
if(Sign_compared<0 && Num>1)
    for i = 1:Num
        %%%% In the last version, if Sign_compared>-1, choose the signal in
        %%%% the T60=0s, else in the select T60
%         Ori_s_temp = audioread(strcat(dataPath,'0/',num2str(Angle(i)),'/Mic',num2str(Sign_compared),'.wav'));
        %%%% In this version, if Sign_compared<0, chosse the original
        %%%% signal, else choose the select T60
        Ori_s_temp = audioread(strcat('simu_data/oral/',num2str(Angle(i)),'.wav'));
        Ori_s(:,i) = Ori_s_temp(1:len);
    end
end
if(SetupStruc.AddNoiseFlag==1)
    av_pow = mean( sum(s.^2,1)/len );       % Average mic power across all received signals.
    sigma_noise = sqrt( av_pow/(10^(SetupStruc.NoiseSNR/10)) );		% st. dev. of white noise component to achieve desired SNR.
    s = s + sigma_noise*randn(size(s));      % Add some random noise
end
% autoPlot(Ori_s,fs,-1);   %%%the 3rd parameter is used to control whether adjusting the range of y axis when equal -1
% autoPlot(s(:,1),fs);
return;