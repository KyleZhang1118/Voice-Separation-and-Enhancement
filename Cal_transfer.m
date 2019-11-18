function [Transfer,SetupStruc] = Cal_transfer(SetupStruc,method)
c = SetupStruc.c;
d = SetupStruc.Size;
fs = SetupStruc.fs;
eval(strcat('K = SetupStruc.',method,'.K;'));
num = SetupStruc.Channel_Num;
Sign = SetupStruc.Sign_Mid;
Angle = SetupStruc.Angle;
if(Sign == 1)
    alpha = 2*pi/(num-1);
else
    alpha = 2*pi/num;
end
Transfer = zeros(360,num,K/2+1);
omiga = zeros(360,num);
for f = 1:K/2+1
    if(Sign == 1)        
        for i = 1:num-1
            delay_d = d*cos(-alpha/2+i*alpha-[0:359]'*pi/180)/c;%the time the ith mic lead to the center mic
            omiga(:,i+1) = 2*pi*fs*(f-1)/K*delay_d;
        end
    else
        for i = 1:num
            delay_d = d*cos(-alpha/2+i*alpha-[0:359]'*pi/180)/c;%the time the ith mic lead to the center mic
            omiga(:,i) = 2*pi*fs*(f-1)/K*delay_d;
        end
    end
    Transfer(:,:,f)= exp(1i*omiga);
end
SetupStruc.RIR = Transfer;
Transfer = permute(Transfer(round(Angle)+1,:,:),[3 2 1]);
SetupStruc.Transfer = Transfer;
return;