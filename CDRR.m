function DRR = CDRR(S,X,fs)
% S为系统输入信号
% X为系统输出信号
M = fs*0.75;         % RIR阶数(时间/s)
miu = 0.001;           % 收敛步长(标量)要求大于0,小于X的相关矩阵最大特征值的倒数 % miu = 1./max(eig(X*X.')) - 0.0001;
%%%%%%%%  计算DRR_m
RIR = RIR_LMS(S,X,M,miu);  
thresh = 0.03;     % 早期混响和晚期混响的分界线（时间 = 30 ms）
% thresh = 0.0025;
idx = round(thresh*fs);
D = sum(RIR(1:idx).^2);
R = sum(RIR(idx+1:end).^2);
DRR = 10*log10(D/R);
% autoPlot(RIR,fs)
return