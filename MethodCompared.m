close all;
clc;clear;
folderPath = '/Results';
dirFol = dir(folderPath);
j_num=0;
Metrics_num = 3;
sign_ignore = 0;
for i = 1:length(dirFol)  %%%% Method circulation
    if( isequal(dirFol(i).name,'.')||...
        isequal(dirFol(i).name,'..')||...
        isequal(dirFol(i).name,'All')||...
        ~dirFol(i).isdir)
        sign_ignore = sign_ignore+1;
        continue;
    end
    method = dirFol(i).name;
    Method{i-sign_ignore} = strrep(method,'_',' ');
    dirMet = dir(fullfile(folderPath,method));
    for j = 1:4%length(dirMet)  %%%%% T60 circulation
        if( isequal(dirMet(j).name,'.')||...
            isequal(dirMet(j).name,'..')||...
            ~dirMet(j).isdir)
            continue;
        end
        j_num = j_num+1;
        subpath = fullfile(folderPath,method,dirMet(j).name,'*.txt');
        dat = dir(subpath);
        SDR = [];SIR = [];SAR = [];
        xAxis = [];
        for k = 1:length(dat)
            dataPath = fullfile(folderPath,method,dirMet(j).name,dat(k).name);
            fileID = fopen(dataPath);
            A = textscan(fileID,'%n%n%n');
            SDR = [SDR;mean(rmmissing(A{1}))];
            SIR = [SIR;mean(rmmissing(A{2}))];
            SAR = [SAR;mean(rmmissing(A{3}))];
            xAxis = [xAxis;str2double(dat(k).name(8:end-4))];
            fclose(fileID);
        end
        [B,I] = sort(xAxis);
        xAxis = xAxis(I);
        SDR = SDR(I);
        SIR = SIR(I);
        SAR = SAR(I);
        if(Metrics_num>=1)
            figure((j-3)*Metrics_num+1)
            plot(xAxis,SDR)
            hold on
        end
        if(Metrics_num>=2)
            figure((j-3)*Metrics_num+2)
            plot(xAxis,SIR)
            hold on
        end
        if(Metrics_num>=3)
            figure((j-3)*Metrics_num+3)
            plot(xAxis,SAR)
            hold on
        end
    end
end
method_num = length(Method);
leS = '';
for i=1:method_num
    if(i>1)
        leS = strcat(leS,',');
    end
    leS = strcat(leS,'''',Method{i},'''');
end
j_num = j_num/method_num*Metrics_num;
for i =1:j_num
    if(mod(i,Metrics_num)==1 || Metrics_num==1)
        figure(i)
        xlabel('Angle/degree')
        ylabel('SDR/dB')
        eval(strcat('legend(',leS,')'));
        title(strcat('Methods in  Room with T60 =',dirMet((i-1)/Metrics_num+3).name(5:end),'s'))
        grid on
        if(Metrics_num==1) continue; end
    elseif(mod(i,Metrics_num)==2 || Metrics_num==2)
        figure(i)
        xlabel('Angle/degree')
        ylabel('SIR/dB')
        eval(strcat('legend(',leS,')'));
        title(strcat('Methods in  Room with T60 =',dirMet((i-2)/Metrics_num+3).name(5:end),'s'))
        grid on
        if(Metrics_num==2) continue;end
    elseif(mod(i,Metrics_num)==0)
        figure(i)
        xlabel('Angle/degree')
        ylabel('SAR/dB')
        eval(strcat('legend(',leS,')'));
        title(strcat('Methods in  Room with T60 =',dirMet(i/Metrics_num+2).name(5:end),'s'))
        grid on
    end
end
        
        
        