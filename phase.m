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
filename1 = '6_sonar_tr4x4_pos_b41.csv';
T1 = readtable(filename1,'VariableNamingRule','preserve');
S1 = table2struct(T1);
filename2 = '1_meas_item_list_rev.xlsx';
T2 = readtable(filename2,'VariableNamingRule','preserve');
T2.Properties.VariableNames{'・データは2回ずつ取得'} = 'Var1';
T2.Properties.VariableNames{'※湿度はtesto温湿度計を使用した参考値'} = 'Var10';
S2 = table2struct(T2);
%%
[file,path] = uigetfile('.\2_tof_data\*.wav','Select a File');
pos_tof=strfind(path,'tof_data');
pos_frame=strfind(path,'frame');
pos_freq=strfind(path,'kHz');
path_name1=path(pos_tof+9:pos_frame-2);
path_name2=path(pos_frame:end-1);
Fc=str2double(path(pos_freq-2:pos_freq-1));% [kHz]
Fc=1000*Fc;% [Hz]

%% 実験No
No=str2double(path_name1(1:2));

%% マイクの位置
pos=str2double(path_name2(8:9));
pos=rem(pos,16);
if pos==0
    pos=16;
end

oppos=str2double(file(3:4));

[wdata,Fs] = audioread([path file]) ;

dnum=length(wdata);
xtime = (0:dnum-1)*(1/Fs)*1000;%[mx]

xmin=0; xmax=15;% 0-15ms
ymin=0; ymax=1.1;

%% 温度
if No<9
    No=1;
elseif 8 < No&&No <17
    No=9;
elseif 16 < No&&No <25
    No=17;
else
    No=25;
end
T=tem((S2(No).Var10+S2(No).Var11)/2);


f1=figure(1);
set(f1,'Position', [500 100 1200 800])
f2=figure(2);
set(f2,'Position', [600 100 1200 800])
f3=figure(3);
set(f3,'Position', [700 100 1200 800])
f4=figure(4);
set(f4,'Position', [800 100 1200 800])
f5=figure(5);
set(f5,'Position', [900 100 1200 800])

%% グラフの切り取り範囲
if 1 <= rem(str2double(path_name1(1:2)),8) && rem(str2double(path_name1(1:2)),8) <=4
    timew=350;
else
    timew=150;
end


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

% 時間窓設定
idxmbase=fix((leng(S1(pos).xpos,S1(pos).ypos,S1(oppos).xpos,S1(oppos).ypos)*192000/T+0.18*192));
if idxmbase < 0
idxmbase = 1;
end
time.s=xtime(fix(idxmbase-0.18*192));
time.e=xtime(fix(idxmbase-0.18*192+timew));

wdmpoint = fix(idxmbase+96);
%% 順方向
for i=1:2
    path(pos_freq-10:pos_freq-9) = num2str(No+i-1,'%02d');
    path_name1(1:2) = num2str(No+i-1,'%02d');
    for j=1:3
        path(end-2:end-1)= num2str(pos+(j-1)*16,'%02d');
        path_name2(8:9)= num2str(pos+(j-1)*16,'%02d');
        file(9:10) = num2str(pos+(j-1)*16,'%02d');
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


        figure(1);
        hold on
        pl = plot(xtime,angle(wdm)); set(pl,ps)
        xlim([time.s time.e]);ylim([-3.99 3.99]);
        ylabel('Phase[rad]')
        xlabel('Time(ms)')
        xline(xtime(wdmpoint),'-r');
        xline(xtime(idxmbase),'-b');
        hold off
        set(gca,ax)

        figure(2);
        hold on
        pl = plot(xtime,abs(wdm)); set(pl,ps)
        xline(xtime(idxmbase),'-r');
        hold off
        xlim([time.s time.e]);ylim([0 ymax]);
        ylabel('Amplitude')
        xlabel('Time(ms)')
        set(gca,ax)

        figure(3);axis equal
        wdmav=0;
        for k=0:9
        wdmav=wdm(wdmpoint+k)+wdmav;
        end
        wdmav=wdmav/10;
        hold on
        pl = plot(wdmav,'-o'); set(pl,ps)
        set(gca,ax)
        xlim([-1 1]);ylim([-1 1])
        plot([-1 1],[0 0],'k')
        plot([0 0],[-1 1],'k')
        hold off
    end
end


%% 逆方向
file(3:4)=num2str(pos,'%02d');
for i=1:2
    path(pos_freq-10:pos_freq-9) = num2str(No+i-1,'%02d');
    path_name1(1:2) = num2str(No+i-1,'%02d');
    for j=1:3
        path(end-2:end-1)= num2str(oppos+(j-1)*16,'%02d');
        path_name2(8:9)= num2str(oppos+(j-1)*16,'%02d');
        file(9:10) = num2str(oppos+(j-1)*16,'%02d');
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


        figure(4);
        hold on
        pl = plot(xtime,angle(wdm)); set(pl,ps)
        xlim([time.s time.e]);ylim([-3.99 3.99]);
        ylabel('Phase[rad]')
        xlabel('Time(ms)')
        xline(xtime(wdmpoint),'-r');
        xline(xtime(idxmbase),'-b');
        hold off
        set(gca,ax)

        figure(5);
        hold on
        pl = plot(xtime,abs(wdm)); set(pl,ps)
        xline(xtime(idxmbase),'-r');
        hold off
        xlim([time.s time.e]);ylim([0 ymax]);
        ylabel('Amplitude')
        xlabel('Time(ms)')
        set(gca,ax)

        figure(3);axis equal
        wdmav=0;
        for k=0:9
        wdmav=wdm(wdmpoint+k)+wdmav;
        end
        wdmav=wdmav/10;
        hold on
        pl = plot(wdmav,'-o','MarkerEdgeColor','red'); set(pl,ps)
        set(gca,ax)
        hold off
    end
end

%% templeture
function T=tem(t)
    T=331.5+0.6*t;
end

%% length
function len=leng(x1,y1,x2,y2)
    l=(x1-x2)^2+(y1-y2)^2;
    len=sqrt(l);
end

