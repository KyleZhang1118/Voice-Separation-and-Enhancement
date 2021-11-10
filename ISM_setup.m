function [SetupStruc] = ISM_setup()
%%%%%%   the environment coefficients
SetupStruc.T60 = 0;      %T60, Options: 0, 0.3, 0.6, 0.9
SetupStruc.c = 340;
SetupStruc.Sign_compared = -1;   % choosing the standard signal, non-negtive corresponding to the signal in the select T60, or the original signal
%%%%%%    the microphone coefficients
SetupStruc.Size = 4.35*0.01;
SetupStruc.Channel_Num = 7;      %the number of mics
SetupStruc.Sign_Mid = 1;         %the exist of the center mic
SetupStruc.AddNoiseFlag = 0;
SetupStruc.NoiseSNR = 30;
SetupStruc.Height = 0; %set 0 use the sameheight data, set 1 use the differentheight data

