clc;clear;
close all;
method = 'LCMV_Search';
folderPath = fullfile('./Results',method);
dirOri = dir(folderPath);
for i = 1:length(dirOri)
    if( isequal(dirOri(i).name,'.')||...
        isequal(dirOri(i).name,'..')||...
        ~dirOri(i).isdir)
        continue;
    end
    subpath = fullfile(folderPath, dirOri(i).name, '*.txt');
    dat = dir(subpath);
    SDR = [];SIR = [];SAR = [];
    xAxis = [];Num = 0;Num_n = 0;
    for j = 1:length(dat)
        dataPath = fullfile(folderPath, dirOri(i).name, dat(j).name);
        fileID = fopen(dataPath);
        A = textscan(fileID,'%n%n%n');
        A_SDR= A{1};
        Num = Num+length(A_SDR);
        sign_NaN = find(isnan(A_SDR)==1);
        if(~isempty(sign_NaN))
            Num_NaN = length(sign_NaN);
            Num_n = Num_n+Num_NaN;
            fprintf('%s in T60 = %ss, angle = %sбу\n',method,dirOri(i).name(5:end),dat(j).name(8:end-4));
            fprintf('Order:');
            fprintf(' %d',sign_NaN(1:2:end));
            fprintf(' fails %d times, %.2f%%.\n',...
                Num_NaN/2,Num_NaN/length(A_SDR)*100);
        end
        SDR = [SDR;mean(rmmissing(A{1}))];
        SIR = [SIR;mean(rmmissing(A{2}))];
        SAR = [SAR;mean(rmmissing(A{3}))];
        xAxis = [xAxis;str2double(dat(j).name(8:end-4))];
        fclose(fileID);
    end
    if(Num_n>0)
        fprintf('%s in T60 = %ss fails %d times, %.2f%%.\n\n',method,dirOri(i).name(5:end),Num_n/2,Num_n/Num*100);
    end
    [B,I] = sort(xAxis);
    xAxis = xAxis(I);
    SDR = SDR(I);
    SIR = SIR(I);
    SAR = SAR(I);
    figure
    plot(xAxis,SDR,'-')
    hold on
    plot(xAxis,SIR,'.-')
    hold on
    plot(xAxis,SAR,'--')
    xlabel('Angle/бу')
    ylabel('SNR/dB')
    legend('SDR','SIR','SAR')
    title(strcat(method,' in  Room with T60 =',dirOri(i).name(5:end),'s'))
    grid on
end
