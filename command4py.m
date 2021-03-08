clc;clear;
close all;
Angle = [90,270];                               %Angle,Options:90,135,180,270 plus:0,45,225,315
if(exist('SetupStruc','var'))
    [s,sOri,Unmix_s,~] = readData(Angle,ISM_setup);
else
    [s,sOri,Unmix_s,SetupStruc] = readData(Angle,ISM_setup);    % 's' is the muli-channel mixture signal,size:Length*7,
end
audiowrite('s.wav',s,16000)
Re.method = {'FastMNMF_temp','FastMNMF'};
Re.FastMNMF_temp.W = 0;
Re.FastMNMF_temp.S = data;
Re.FastMNMF.W = 0;
Re.FastMNMF.S = data;
[Me] = Cal_metrics(Re,Unmix_s,sOri,SetupStruc);