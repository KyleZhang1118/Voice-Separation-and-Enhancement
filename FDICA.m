function [Y,W,A] = FDICA(X,W,A,F_i)

%First step: Initialization
% step_size=1e-5;                   %the step distance
step_size = 0.1;
% C_r=0.0002;                             % correlation
max_iteration=1000;                   % the max number of iteration

% if(F_i==33) 
%     step_size = step_size/10;
% end
%%%%%%%%%%
frame_N = size(X,2);
% if(F_i==16)
%     a =1;
% end
sign_size = 0;
y_f = W*X;
norm = max(abs(y_f),[],2);
% if(norm>10)
%     norm = repmat(norm,1,size(W,2));
%     W = W./norm;
%     y_f = W*X;
% end

for n_i=1:max_iteration
%         y_fun = power(1+exp(-y_r),-1)+j*power(1+exp(-y_i),-1);
%         W_o1 = step_size*(diag(diag(y_fun*y_f'/frame_N))-(y_fun*y_f'/frame_N))*W_o{K_i}+W_o{K_i};
    y_f = W*X;
    y_r = real(y_f);
    y_i = imag(y_f);
    y_fun = tanh(y_r)+1i*tanh(y_i);
    sign = abs(det(eye(size(X,1))-(y_fun*y_f'/frame_N)));
    A(n_i,F_i-1) = sign;
    if(sign>100 && sign_size ==0)
        step_size = step_size/10;
        sign_size =1;
    end
    if(sign>1000 && sign_size ==1)
        step_size = step_size/10;
        sign_size =2;
    end
    if(sign>10000 && sign_size ==2)
        step_size = step_size/10;
        sign_size =3;
    end
    W1 = step_size*(eye(size(X,1))-(y_fun*y_f'/frame_N))*W+W;
    W = W1;

    if(sign<1e-4)
        break;
    end
%     A(n_i,F_i-1) = abs(det(eye(size(X,1))-(y_fun*y_f'/frame_N)));
%     A(n_i,K_i) = abs(det(eye(size(s,2))-(y_fun*y_f'/frame_N)));
%     if(abs(A(n_i,K_i))<0.001) break; end;
end
A(end,F_i-1) = n_i;
Y = y_f;
% Y = Project_back(W,y_f);
return;