function [Y,W,SetupStruc] = Process_ILRMA_woDR(s,Transfer,SetupStruc)
K = SetupStruc.ILRMA_woDR.K;
hop = SetupStruc.ILRMA_woDR.hop;
win = hanning(K,'periodic');
win = win/sqrt(sum(win(1:hop:K).^2));
SetupStruc.ILRMA_woDR.win = win;  % Preserve 'win' in 'SetupStruc'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = size(s,2);
for i = 1:N
    X(:,:,i) = fft(enframe(s(:,i),win,hop)');
end
frame_N = size(X,2);
K_m = K/2+1;
Num = size(Transfer,3);
Y = zeros((frame_N-1)*hop+K,N);
Y_f = zeros(N,frame_N,K);
%%%%%%%%%%%%%%%%%%%%%%%%%% Obtain processing matrix 'W'
X_sp = zeros(N,frame_N,K_m);
W_ILRMA_woDR = zeros(N,N,K_m);
V_sp = zeros(N,N,K_m);
theta = 10^-6;
for i = 1:K_m
    X_f = permute(X(i,:,:),[3 2 1]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Initialize W by PCA
%     [E,D] = PCA(X_f,1,Num);
%     V = sqrt(D)\E';
    V = eye(N);
    V_sp(:,:,i) = V;
    X_sp(:,:,i) = V*X_f;
    %%%%%%%%%%%% Adjust amplitude of 'w'
    W_o = eye(N);
    y_f = W_o*V*X_f;
    W_ILRMA_woDR(:,:,i) = W_o;
    Y_f(:,:,i) = y_f;    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% ILRMA iterations
epsi = 1e-6;
max_iteration = 200;
L = 2;  %%%%% the number of NMF basis
T = max(rand(K_m,L,N),epsi);
V = max(rand(L,frame_N,N),epsi);
R = zeros(K_m,frame_N,N);
for i = 1:N
    R(:,:,i) = T(:,:,i)*V(:,:,i);
end
P = permute(Y_f(:,:,1:K_m),[3,2,1]);
P = abs(P).^2;
% K_sp = sqrt(K_m);
pObj = inf;
% step_size = 0.01;
A = zeros(1001,2)-1; %%%% Show the decrease of the value of cost funtion, ILRMA max iterations 1000
for iteration = 1:max_iteration
    %%%%% NMF
    for i = 1:N
        T(:,:,i) = max(T(:,:,i).*sqrt(P(:,:,i).*(R(:,:,i).^(-2))*V(:,:,i)'./((R(:,:,i).^(-1))*V(:,:,i)')),epsi);
        R(:,:,i) = T(:,:,i)*V(:,:,i);
        V(:,:,i) = max(V(:,:,i).*sqrt(T(:,:,i)'*(P(:,:,i).*(R(:,:,i).^(-2)))./(T(:,:,i)'*(R(:,:,i).^(-1)))),epsi);
        R(:,:,i) = T(:,:,i)*V(:,:,i);
    end
    %%%%% IVA
    dlw = 0;
    for i = 1:K_m
        W = W_ILRMA_woDR(:,:,i);
        X_f = X_sp(:,:,i);
        dlw = dlw +log(abs(det(W))+epsi); 
        %%%%% IVA  Problem: IVA cant work well because the value of step_size is sensitive when work with NMF, without
        %%%%% dimensionality reduction and whiten
%         y_f = Y_f(:,:,i);
%         G_ = permute(R(i,:,:),[3 2 1]);
%         y_fun = y_f./(G_+epsi);
%         core = eye(size(W))-y_fun*y_f'/frame_N;
%         W = W+step_size*core*W;
        %%%%%  AuxIVA
        for i_n = 1:N
%             G_ = Y_k(i_n,:).^-1;
            G_ = permute(R(i,:,i_n),[3 2 1]);
            G_ = repmat(G_+epsi,N,1);
            Vk = (X_f./G_)*X_f'/frame_N;
            if rcond(Vk)<theta
                Vk = Vk+eye(N)*max(eig(Vk))*theta;
            end
            wk = inv(W*Vk);
            wk = wk(:,i_n);
            wk = wk/(sqrt(wk'*Vk*wk)+epsi);
            W(i_n,:) = wk';
        end
        %%%%%
        W_ILRMA_woDR(:,:,i) = W;
        Y_f(:,:,i) = W*X_f;
    end
    P = permute(Y_f(:,:,1:K_m),[3,2,1]);
    P = abs(P).^2;
    Obj = ((sum(sum(sum(log(R+epsi))))+sum(sum(sum(P./R))))/frame_N-2*dlw)/(N*K_m);
    dObj = pObj-Obj;
    pObj = Obj;
    A(iteration,:) = [Obj,abs(dObj)/abs(Obj)];
    if(abs(dObj)/abs(Obj)<theta)
       break;
    end
    for i = 1:N
        lamda = sqrt(sum(sum(P(:,:,i)))/(frame_N*K_m));
        W_ILRMA_woDR(i,:,:) = W_ILRMA_woDR(i,:,:)/lamda;
        P(:,:,i) = P(:,:,i)/(lamda^2+epsi);
        R(:,:,i) = R(:,:,i)/(lamda^2+epsi);
        T(:,:,i) = T(:,:,i)/(lamda^2+epsi);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i = 1:Num
%     figure
%     imagesc(flipud(R(:,:,i)/max(max(R(:,:,i))+epsi)))    
% end
% for i = 1:Num
%     figure
%     imagesc(flipud(P(:,:,i)/max(max(P(:,:,i))+epsi)))
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Post processing
W = zeros(N,N,K_m);
Y_f(:,:,1) = zeros(N,frame_N);
for i = 2:K_m
    W_inv = pinv(W_ILRMA_woDR(:,:,i)*V_sp(:,:,i));
    for ii = 1:N
        Y_f(ii,:,i) = Y_f(ii,:,i)*W_inv(1,ii);
        W_ILRMA_woDR(ii,:,i) = W_ILRMA_woDR(ii,:,i)*W_inv(1,ii);
    end
    W(:,:,i) = W_ILRMA_woDR(:,:,i)*V_sp(:,:,i);     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(i~=K_m)
        Y_f(:,:,K+2-i) = conj(Y_f(:,:,i));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Recover signals
if(K/hop==2)
    win = ones(K,1);
end
for i = 1:N
    y_temp = permute(Y_f(i,:,:),[3 2 1]);
    Y(:,i) = overlapadd(real(ifft(y_temp))',win,hop);
end
% autoPlot(Y,16000)
Power = sum(Y.^2);
[~,Order] = sort(Power,'descend');
Y = Y(:,Order(1:Num));
W = W(Order(1:Num),:,:);
return;