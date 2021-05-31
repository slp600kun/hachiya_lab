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
% スプレッドシート
filename = '6_sonar_tr4x4_pos_b41.csv';
T = readtable(filename);
S = table2struct(T);
%%
[file,path] = uigetfile('.\2_tof_data\*.wav','Select a File');
pos_tof=strfind(path,'tof_data');
pos_frame=strfind(path,'frame');
pos_freq=strfind(path,'kHz');
path_name1=path(pos_tof+9:pos_frame-2);
path_name2=path(pos_frame:end-1);
Fc=str2double(path(pos_freq-2:pos_freq-1));% [kHz]
Fc=1000*Fc;% [Hz]
%% マイクの位置
pos=str2double(path(pos_frame+7:pos_frame+8));
pos=rem(str2double(path(pos_frame+7:pos_frame+8)),16);
if pos==0
    pos=16;
end

%%
ch_str = num2str(1,'%02d') ;
file(3:4) = ch_str ;
pos_str = num2str(pos,'%02d');
[wdata,Fs] = audioread([path file]) ;

dnum=length(wdata);
xtime = (0:dnum-1)*(1/Fs)*1000;%[mx]

ymin=0; ymax=1.5;
t1=zeros(1,16);
t2=zeros(1,16);

f1=figure(1);
set(f1,'Position', [700 500 1200 800])

f2=figure(2);
set(f2,'Position', [600 300 1200 800])

f3=figure(3);
set(f3,'Position', [500 100 1200 800])

f4=figure(4);
set(f4,'Position', [700 500 1200 800])

f5=figure(5);
set(f5,'Position', [600 300 1200 800])

f6=figure(6);
set(f6,'Position', [500 100 1200 800])

%% Filter 設定
iftr=1;

if iftr==1
    Fstop = 7000;
    Fpass = 8000;
    Astop = 30;
    Apass = 0.5;
    
    d = designfilt('highpassfir','StopbandFrequency',Fstop, ...
        'PassbandFrequency',Fpass,'StopbandAttenuation',Astop, ...
        'PassbandRipple',Apass,'SampleRate',Fs,'DesignMethod','equiripple');
    
    delay = round(mean(grpdelay(d)))+1;
end

%% tpを求める
for ix=1:16
    ch_str = num2str(ix,'%02d');
    if ix ~= 1
        file(3:4) = ch_str ;
    end
    
    % Filter 処理
    if iftr==1
        [wdata_tmp,Fs] = audioread([path file]) ;
        wdata=filter(d,wdata_tmp) ;
        wdata(1:end-delay+1) = wdata(delay:end);
    else
        [wdata,Fs] = audioread([path file]);
    end
    
    [b, a] = demod(wdata, Fc, Fs, 'qam');
    wdm = complex(a, b);
    
    % 伝播時間
    [md,mpos]=max(abs(wdm));
    minpos=find(abs(wdm)>(md*0.15), 1 );
    t1(1,ix)=xtime(1,minpos);
    T1=t1(1,ix);
    
    figure(1);subplot(4,4,ix)
    pl = plot(xtime,wdata); set(pl,ps)
    xlim([T1,T1+1.1]);ylim([-ymax ymax]);
    ylabel('Amplitude')
    
    if ix==2
        tp=title(path_name1);set(tp,tx)
    elseif ix==3
            tp=title(path_name2);set(tp,tx)
    elseif ix >= 13
            xlabel('Time(ms)')
    end
    text(T1+0.85,1.4,['CH : ' ch_str],'FontSize',9)
    set(gca,ax)xtim

    figure(2);subplot(4,4,ix)
    pl = plot(xtime,abs(wdm)); set(pl,ps)
    xlim([T1,T1+1.1]);ylim([0 1.2]);
    ylabel('Amplitude')
    if ix==2
        tp=title(path_name1);set(tp,tx)
    elseif ix==3
            tp=title(path_name2);set(tp,tx)
    elseif ix >= 13
            xlabel('Time(ms)')
    end
    text(T1+0.85,1.1,['CH : ' ch_str],'FontSize',9)
    set(gca,ax)
    
    figure(3);subplot(4,4,ix)
    pl = plot(xtime,angle(wdm)); set(pl,ps)
    xlim([T1,T1+1.1]);ylim([-3.99 3.99]);
    ylabel('Phase[rad]')
    if ix==2
        tp=title(path_name1);set(tp,tx)
    elseif ix==3
            tp=title(path_name2);set(tp,tx)
    elseif ix >= 13
            xlabel('Time(ms)')
    end
    text(T1+0.85,3.5,['CH : ' ch_str],'FontSize',9)
    set(gca,ax)
    
end


for iy=1:16
    ch_str = num2str(iy,'%02d');
    path(end-2:end-1)= ch_str;
    file(3:4)=pos_str;
    file(9:10)=ch_str;
    [wdata,Fs] = audioread([path file]) ;
    dnum=length(wdata);
    xtime = (0:dnum-1)*(1/Fs)*1000;%[mx]    
    % Filter 処理
    if iftr==1
        [wdata_tmp,Fs] = audioread([path file]) ;
        wdata=filter(d,wdata_tmp) ;
        wdata(1:end-delay+1) = wdata(delay:end);
    else
        [wdata,Fs] = audioread([path file]);
    end
    
    [b, a] = demod(wdata, Fc, Fs, 'qam');
    wdm = complex(a, b);
    % 伝播時間
    [md,mpos]=max(abs(wdm));
    minpos=find(abs(wdm)>(md*0.15), 1 );
    t2(1,iy)=xtime(1,minpos);
    T2=t2(1,iy);
    
    
    figure(4);subplot(4,4,iy)
    pl = plot(xtime,wdata); set(pl,ps)
    xlim([T2,T2+1.1]);ylim([-ymax ymax]);
    ylabel('Amplitude')
    
    if ix==2
        tp=title(path_name1);set(tp,tx)
    elseif ix >= 13
            xlabel('Time(ms)')
    end
    text(T2+0.85,1.4,['CH : ' ch_str],'FontSize',9)
    set(gca,ax)

    figure(5);subplot(4,4,iy)
    pl = plot(xtime,abs(wdm)); set(pl,ps)
    xlim([T2,T2+1.1]);ylim([0 1.2]);
    ylabel('Amplitude')
    if ix==2
        tp=title(path_name1);set(tp,tx)
    elseif ix >= 13
            xlabel('Time(ms)')
    end
    text(T2+0.85,1.1,['CH : ' ch_str],'FontSize',9)
    set(gca,ax)
    
    figure(6);subplot(4,4,iy)
    pl = plot(xtime,angle(wdm)); set(pl,ps)
    xlim([T2,T2+1.1]);ylim([-3.99 3.99]);
    ylabel('Phase[rad]')
    if ix==2
        tp=title(path_name1);set(tp,tx)
    elseif ix >= 13
            xlabel('Time(ms)')
    end
    text(T2+0.85,3.5,['CH : ' ch_str],'FontSize',9)
    set(gca,ax)
 
end

l1=zeros(1,16);
for i=1:16
    l1(1,i)=leng(S(pos).xpos,S(pos).ypos,S(i).xpos,S(i).ypos);
end

c1=zeros(1,16);
for iz=1:16
    c1(1,iz)=sound(t1(1,iz),t2(1,iz ),l1(1,iz));
end

%% length
function len=leng(x1,y1,x2,y2)
    l=(x1-x2).^2+(y1-y2).^2;
    len=sqrt(l);
end

%% speed of sound
function c=sound(tP,tM,l)
    c=2*l/(tP/1000+tM/1000);
end