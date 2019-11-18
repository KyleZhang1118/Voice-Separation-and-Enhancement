function [] = autoWrite(S,fs,method)
% 'S'£ºLen * Num
if(~exist('method','var'))
    method = [];
end
num = size(S,2);
S = S/max(max(abs(S)));
for i = 1:num
    eval(strcat('audiowrite(''S_',method,num2str(i),'.wav''',',S(:,',num2str(i),'),','fs);'));
end

return