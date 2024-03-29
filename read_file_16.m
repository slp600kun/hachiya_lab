%
% wave ファイルのプロット
%
%%
clearvars
close all

% 図の基本設定
ax.LineWidth=1.2;ax.FontSize=12;
ps.Color='black';ps.LineWidth=1.5;
tx.Interpreter='none';

%%
[file,path] = uigetfile('.\2_tof_data\*.wav','Select a File');
pos_tof=strfind(path,'tof_data');
pos_frame=strfind(path,'frame');
pos_freq=strfind(path,'kHz');
path_name1=path(pos_tof+9:pos_frame-2);
path_name2=path(pos_frame:end-1);
Fc=str2num(path(pos_freq-2:pos_freq-1));% [kHz]
Fc=1000*Fc;% [Hz]

file(3:4) = num2str(1,'%02d');
[wdata,Fs] = audioread([path file]);

dnum=length(wdata);
xtime = (0:dnum-1)*(1/Fs)*1000;%[mx]

xmin=0; xmax=15;% 0-15ms
ymin=0; ymax=1;

%%
for ix=1:16
    if ix ~= 1
        file(3:4) = num2str(ix,'%02d');
        [wdata,Fs] = audioread([path file]);
    end
    
    [b, a] = demod(wdata, Fc, Fs, 'qam');
    wdm = complex(a, b);

    figure
    %figure(ix)
    subplot(3,1,1)
    pl = plot(xtime,wdata); set(pl,ps)
    xlim([xmin xmax]);ylim([-ymax ymax]);
    ylabel('Amplitude')
    tp=title({path_name1; path_name2; file});set(tp,tx)
    set(gca,ax)
    
    subplot(3,1,2)
    pl = plot(xtime,abs(wdm)); set(pl,ps)
    xlim([xmin xmax]);ylim([0 ymax]);
    ylabel('Amplitude')
    set(gca,ax)
    
    subplot(3,1,3)
    pl = plot(xtime,angle(wdm)); set(pl,ps)
    xlim([xmin xmax]);ylim([-3.2 3.2]);
    xlabel('Time[ms]');ylabel('phase(rad)')
    set(gca,ax)
    
end