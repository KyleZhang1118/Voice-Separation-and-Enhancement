function [Y,SetupStruc] = SimuOnline(s,Transfer,SetupStruc,method)
%%%%%模拟在线对原信号进行分帧并构成块，每一帧信号都进行加窗，输入算法核
%%%%%将核返回的块进行去窗拼成时域信号，拼接上一块
%%%%%
%%%%%算法核定义，接收已加窗的语音帧构成的块，针对该块进行计算返回该块
K = SetupStruc.K;
hop = SetupStruc.hop;
L = SetupStruc.L;
win = hanning(K,'periodic');
win = win/sqrt(sum(win(1:hop:K).^2));
SetupStruc.win = win;  %win窗反传保留
N = size(s,2);
Num = size(Transfer,3);
Y = zeros(size(s,1),Num);
Block = zeros(K,N,L); %%%% 帧长x通道x帧数
a = 1;
aBlock = 1;
l = 1;
if(isfield(SetupStruc,'InL'))    %%%%判断该函数是否需要进行第一次长计算以初始化参数
    Block1 = zeros(K,N,SetupStruc.InL);
    if(isfield(SetupStruc,'unS'))
        SetupStruc.Block = zeros(K,N,SetupStruc.InL,Num);
    end
    for i = 1:SetupStruc.InL
        Block1(:,:,i) = s(a:a+K-1,:);
        if(isfield(SetupStruc,'unS'))
            for j = 1:Num
                SetupStruc.Block(:,:,i,j) = SetupStruc.unS(a:a+K-1,:,j);
            end
        end
        a = a+hop;
        if(a+K-1>size(s,1))
            break;
        end
    end
    LengthI = (i-1)*hop+K;
    eval(strcat('[y_temp,SetupStruc] = On',method,'(Block1(:,:,1:i),Transfer,SetupStruc);'));
    Y(aBlock:aBlock+LengthI-1,:) = Y(aBlock:aBlock+LengthI-1,:)+y_temp;
    aBlock = a;
end
if(isfield(SetupStruc,'unS'))
    SetupStruc.Block = zeros(K,N,L,Num);
end
while(a+K<size(s,1))
    Block(:,:,l) = s(a:a+K-1,:); %%%填充块
    if(isfield(SetupStruc,'unS'))
        for j = 1:Num
            SetupStruc.Block(:,:,l,j) = SetupStruc.unS(a:a+K-1,:,j);
        end
    end
    a = a+hop;
    Length = (l-1)*hop+K;
    if(l==L || a+K-1>size(s,1))
        eval(strcat('[y_temp,SetupStruc] = On',method,'(Block(:,:,1:l),Transfer,SetupStruc);'));
        Y(aBlock:aBlock+Length-1,:) = Y(aBlock:aBlock+Length-1,:)+y_temp;
        l = 0;
        aBlock = a;
    end
    l = l+1;
end

end