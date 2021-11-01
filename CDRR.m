function DRR = CDRR(S,X,fs)
% SΪϵͳ�����ź�
% XΪϵͳ����ź�
M = fs*0.75;         % RIR����(ʱ��/s)
miu = 0.001;           % ��������(����)Ҫ�����0,С��X����ؾ����������ֵ�ĵ��� % miu = 1./max(eig(X*X.')) - 0.0001;
%%%%%%%%  ����DRR_m
RIR = RIR_LMS(S,X,M,miu);  
thresh = 0.03;     % ���ڻ�������ڻ���ķֽ��ߣ�ʱ�� = 30 ms��
% thresh = 0.0025;
idx = round(thresh*fs);
D = sum(RIR(1:idx).^2);
R = sum(RIR(idx+1:end).^2);
DRR = 10*log10(D/R);
% autoPlot(RIR,fs)
return