function [Y,SetupStruc] = SimuOnline(s,Transfer,SetupStruc,method)
%%%%% Simulate of online framing and blocking signals,  adding window funcion, putting into the algorithm core
%%%%% Convert the block back from algo core to time domian, splice the last block
%%%%% 
%%%%% The definition of algo core: recieving the block of speech frames, process the block then feedback
K = SetupStruc.K;
hop = SetupStruc.hop;
L = SetupStruc.L;
win = hanning(K,'periodic');
win = win/sqrt(sum(win(1:hop:K).^2));
SetupStruc.win = win;  % Save the window function
N = size(s,2);
Num = size(Transfer,3);
Y = zeros(size(s,1),Num);
Block = zeros(K,N,L); %%%% frame length* channel* frame number
a = 1;
aBlock = 1;
l = 1;
if(isfield(SetupStruc,'InL'))    %%%% Judge if this method needs a long range cal to initialize
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
    Block(:,:,l) = s(a:a+K-1,:); %%% Padding the bolck
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