function [Y,W,SetupStruc] = Process_MVDR_Search(s,Transfer,SetupStruc)
K = SetupStruc.MVDR_Search.K;
hop = SetupStruc.MVDR_Search.hop;
win = hanning(K,'periodic');
win = win/sqrt(sum(win(1:hop:K).^2));
SetupStruc.MVDR_Search.win = win;  % Preserve 'win' in 'SetupStruc'
% K = SetupStruc.MVDR.K;
% hop = SetupStruc.MVDR.hop;
% win = hanning(K,'periodic');
% win = win/sqrt(sum(win(1:hop:K).^2));
% SetupStruc.MVDR.win = win;  % Preserve 'win' in 'SetupStruc'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = size(s,2);
for i = 1:N
    X(:,:,i) = fft(enframe(s(:,i),win,hop)');
end
frame_N = size(X,2);
K_m = K/2+1;
Num = size(Transfer,3);
Y = zeros((frame_N-1)*hop+K,Num);
Y_f = zeros(size(X,1),size(X,2),Num);
WNG = zeros(K_m,Num);
%%%%%%%%%%%%%%%%%%%%%%%%%% Obtain processing matrix 'W'
theta = 10^-10;
WNG_set = 0;
W = zeros(Num,N,K_m);
for i = 2:K_m
    X_f = permute(X(i,:,:),[3 2 1]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% MVDR processing
    W_f = zeros(N,Num);
    Steer = permute(Transfer(i,:,:),[2 3 1]);
    for j = 1:Num
        R = X_f*X_f'/frame_N;
        if rcond(R)<theta
            R = R+eye(N)*min([min(diag(R)) theta]);
        end
        h = Steer(:,j);
        w = R\h/(h'/R*h);
        wng = -10*log10(diag(w'*w));
        e = 0.05;
        if(WNG_set-wng>0)
            sign = 1;
        else 
            sign = 0;
        end
        signM = 1;
        while(sign==1 && abs(wng-WNG_set)>0.5)
            if(wng>WNG_set)
                R = R-e*eye(N);
                e = e/10;
                signM = 0;
            else
                if(signM==1)
                    e = e*2;
                end
                R = R+e*eye(N);
            end
            w = R\h/(h'/R*h);
            wng = -10*log10(diag(w'*w));
        end
%         while(sign==0 && abs(wng-WNG_set)>0.5)
%             if(wng<WNG_set)
%                 R = R+e*eye(N);
%                 e = e/10;
%                 signM = 0;
%             else
%                 if(signM==1)
%                     e = e*2;
%                 end
%                 R = R-e*eye(N);
%             end
%             w = R\h/(h'/R*h);
%             wng = -10*log10(diag(w'*w));
%         end
        W_f(:,j) = w;
    end
    W_f = W_f';
    W(:,:,i) = W_f;   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    WNG(i,:) = -10*log10(diag(W_f*W_f')');
    Y_ = W_f*X_f;
    Y_f(i,:,:) = Y_.';
    if(i~=K_m)
        Y_f(K+2-i,:,:) = Y_';
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Recover signals
if(K/hop==2)
    win = ones(K,1);
end
for i = 1:Num
    Y(:,i) = overlapadd(real(ifft(Y_f(:,:,i)))',win,hop);
end
% autoPlot(WNG,'MVDR_Search');
return;