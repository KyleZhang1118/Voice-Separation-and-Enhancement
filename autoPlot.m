function [] = autoPlot(S,method,fs,para)
% 'S':Len * Num, 'method' is the name of method,if the 2rd para is not char, then shift
if(~ischar(method))
    if(~exist('fs','var'))
        fs = method;
    else 
        para = fs;
        fs = method;
    end
    clear method;
end
num = size(S,2);
po_x = size(S,1);
if(exist('fs','var'))
    x = [1:size(S,1)]/fs;
    po_x = max(x);
end
figure
for i = 1:num
    eval(strcat('subplot(',num2str(num),'1',num2str(i),')'));
    if(exist('x','var'))
        plot(x,S(:,i))
    else
        plot(S(:,i))
    end
    ylim([min(min(S)) max(max(S))])
    xlabel('Time/s')
    if(i==1 && exist('method','var'))
        title(strrep(method,'_',' '))
    end
    if(exist('para','var'))
        if(size(para,2)==3)
            text(po_x,max(max(S))*0.9,num2str(para(i,1)))
            text(po_x,max(max(S))*0.7,num2str(para(i,2)))
            text(po_x,max(max(S))*0.5,num2str(para(i,3)))
            continue;
        end
        text(po_x*0.9,max(max(S))*0.9,num2str(para(i,1)))
        if(size(para,2)>4)
            text(po_x*0.9,max(max(S))*0.7,num2str(para(i,2)))
        end
        if(size(para,2)>3)
            text(po_x,max(max(S))*0.9,num2str(para(i,3)))
            text(po_x,max(max(S))*0.7,num2str(para(i,4)))
            text(po_x,max(max(S))*0.5,num2str(para(i,5)))
        end
    end    
end

return