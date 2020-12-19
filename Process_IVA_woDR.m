function [Y,W,SetupStruc] = Process_IVA_woDR(s,Transfer,SetupStruc)
K = SetupStruc.IVA_woDR.K;
hop = SetupStruc.IVA_woDR.hop;
win = hanning(K,'periodic');
win = win/sqrt(sum(win(1:hop:K).^2));
SetupStruc.IVA_woDR.win = win;  % Preserve 'win' in 'SetupStruc'
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
W_IVA = zeros(N,N,K_m);
V_sp = zeros(N,N,K_m);
theta = 10^-6;
for i = 1:K_m
    X_f = permute(X(i,:,:),[3 2 1]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Initialize W by PCA
%     [E,D] = PCA(X_f,1,N);
%     V = sqrt(D)\E';
    V = eye(N);
    V_sp(:,:,i) = V;
    X_sp(:,:,i) = V*X_f;
    %%%%%%%%%%%% Adjust amplitude of 'w'
    W_o = eye(N);
    y_f = W_o*V*X_f;
    W_IVA(:,:,i) = W_o;
    Y_f(:,:,i) = y_f;    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% IVA iterations
epsi = 1e-6;
step_size = 0.1;
max_iteration = 200;
Y_k = zeros(N,frame_N);
% K_sp = sqrt(K_m);
pObj = inf;
A = zeros(1001,2)-1; %%%% Show the decrease of non-linear correlation, IVA max iterations 1000
for iteration = 1:max_iteration
    for i = 1:N
        y_temp = permute(Y_f(i,:,:),[3 2 1]);
        Y_k(i,:) = sqrt(sum(abs(y_temp(1:K_m,:)).^2))+epsi;
    end
    dlw = 0;
    for i = 1:K_m
        W = W_IVA(:,:,i);
        X_f = X_sp(:,:,i);
        y_f = Y_f(:,:,i);
        y_fun = y_f./Y_k;
%         y_fun = K_sp*tanh(K_sp*Y_k).*y_f./Y_k;
        core = eye(size(W))-y_fun*y_f'/frame_N;
        dlw = dlw +log(abs(det(W))+epsi);
        W = W+step_size*core*W;
        W_IVA(:,:,i) = W;
        Y_f(:,:,i) = W*X_f;
    end
    Obj = (sum(sum(Y_k))/frame_N-2*dlw)/(N*K_m);
    dObj = pObj-Obj;
    pObj = Obj;
    A(iteration,:) = [Obj,abs(dObj)/abs(Obj)];
    if(abs(dObj)/abs(Obj)<theta)
        break;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Post processing
W = zeros(N,N,K_m);
Y_f(:,:,1) = zeros(N,frame_N);
for i = 2:K_m
    W_inv = pinv(W_IVA(:,:,i)*V_sp(:,:,i));
    for ii = 1:N
        Y_f(ii,:,i) = Y_f(ii,:,i)*W_inv(1,ii);
        W_IVA(ii,:,i) = W_IVA(ii,:,i)*W_inv(1,ii);
    end
    W(:,:,i) = W_IVA(:,:,i)*V_sp(:,:,i);     
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