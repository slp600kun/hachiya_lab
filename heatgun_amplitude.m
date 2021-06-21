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

%% ファイル名
[file,path] = uigetfile('.\2_tof_data\*.wav','Select a File');
pos_tof=strfind(path,'tof_data');
pos_frame=strfind(path,'frame');
pos_freq=strfind(path,'kHz');
pos_wn=strfind(path,'wn');
path_name1=path(pos_tof+9:pos_frame-2);
path_name2=path(pos_frame:end-1);

%% 周波数
Fc=str2double(path(pos_freq-2:pos_freq-1));% [kHz]
Fc=1000*Fc;% [Hz]

%% 実験No
No=str2double(path_name1(1:2));

%% マイクの位置
posA = str2double(path_name2(8:9));
posA = rem(posA,16);
if posA == 0
    posA = 16;
end
micpos = [12;11;10;9;16;15;14;13];
posB = micpos(posA);

%% サンプリング周波数
[wdata,Fs] = audioread([path file]) ;
dnum=length(wdata);
xtime = (0:dnum-1)*(1/Fs)*1000;%[mx]
ymin=0; ymax=1.5;

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
T=tem1((S2(No).Var10+S2(No).Var11)/2);

%% グラフデータ
f1=figure(1);
set(f1,'Position', [500 100 1600 800])

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

Amp=zeros(2,3);


for i=1:2
    path(pos_freq-10:pos_freq-9) = num2str(No+i-1,'%02d');
    path_name1(1:2) = num2str(No+i-1,'%02d');
    for j=1:3
        path(end-2:end-1)= num2str(posA+(j-1)*16,'%02d');
        path_name2(8:9)= num2str(posA+(j-1)*16,'%02d');
        file(3:4) = num2str(posB,'%02d');
        file(9:10) = num2str(posA+(j-1)*16,'%02d');
        
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


        % 時間窓設定
        idxmbase=fix((1.5*192000/T-50+0.18*192))+1;
        if idxmbase < 0
            idxmbase = 1;
        end
        time.sppos(i,j)=idxmbase;
        time.sp(i,j)=xtime(idxmbase);
        time.ep(i,j)=xtime(idxmbase+timew);

        % 振幅最大値
        [tymax,tpos]=max(abs(wdm(idxmbase:idxmbase+timew)));
        tpos=idxmbase+tpos;
        tmax.tppos(i,j)=tpos;
        tmax.tpy(i,j)=tymax;
        tmax.tp(i,j)=xtime(tpos);

        % 伝搬時間
        pt.tppos(i,j)=find(abs(wdm(idxmbase:idxmbase+timew))>(tymax*0.1), 1 )+idxmbase;
        pt.tp(i,j)=xtime(pt.tppos(i,j))-0.18;
        
        timesub.sppos(i,j)=pt.tppos(i,j)+fix(0.5*192);
        timesub.sp(i,j)=xtime(timesub.sppos(i,j))-0.18;
        timesub.eppos(i,j)=timesub.sppos(i,j)+fix(0.1*192);
        
        amp=0;
        for k=timesub.sppos(i,j):timesub.eppos(i,j)
            amp=amp+abs(wdm(k));
        end
        amp=amp/(timesub.eppos(i,j)-timesub.sppos(i,j)+1);
        Amp(i,j)=amp;
        
        
        figure(1);subplot(1,2,2)
        hold on
        plot(pt.tp(i,j),Amp(i,j),'o','MarkerEdgeColor','red')
        hold off
        xlim([4.4 4.45]);ylim([0 1.2]);
        xlabel('Time(ms)')
        ylabel('Amplitude')
        grid on
        tp=title({strcat(num2str(No),path_name1(3:end));strcat(num2str(No+1),path_name1(3:end))});set(tp,tx)
        
        %同様の条件の場所
        path(end-2:end-1)= num2str(8-posA+1+(j-1)*16,'%02d');
        path_name2(8:9)= num2str(8-posA+1+(j-1)*16,'%02d');
        file(3:4) = num2str(micpos(8-posA+1),'%02d');
        file(9:10) = num2str(8-posA+1+(j-1)*16,'%02d');
        
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


        % 時間窓設定
        idxmbase=fix((1.5*192000/T-50+0.18*192))+1;
        if idxmbase < 0
            idxmbase = 1;
        end
        time.sppos(i,j)=idxmbase;
        time.sp(i,j)=xtime(idxmbase);
        time.ep(i,j)=xtime(idxmbase+timew);

        % 振幅最大値
        [tymax,tpos]=max(abs(wdm(idxmbase:idxmbase+timew)));
        tpos=idxmbase+tpos;
        tmax.tppos(i,j)=tpos;
        tmax.tpy(i,j)=tymax;
        tmax.tp(i,j)=xtime(tpos);

        % 伝搬時間
        pt.tppos(i,j)=find(abs(wdm(idxmbase:idxmbase+timew))>(tymax*0.1), 1 )+idxmbase;
        pt.tp(i,j)=xtime(pt.tppos(i,j))-0.18;
        
        timesub.sppos(i,j)=pt.tppos(i,j)+fix(0.5*192);
        timesub.sp(i,j)=xtime(timesub.sppos(i,j))-0.18;
        timesub.eppos(i,j)=timesub.sppos(i,j)+fix(0.1*192);
        
        amp=0;
        for k=timesub.sppos(i,j):timesub.eppos(i,j)
            amp=amp+abs(wdm(k));
        end
        amp=amp/(timesub.eppos(i,j)-timesub.sppos(i,j)+1);
        Amp(i,j)=amp;
        
        figure(1);subplot(1,2,2)
        hold on
        plot(pt.tp(i,j),Amp(i,j),'o','MarkerEdgeColor','black')
        hold off
    
    end
end

path = erase(path,path(pos_wn+4:pos_frame-2));
subpos_wn = strfind(path_name1,'wn');
path_name1 = erase(path_name1,path_name1(subpos_wn+4:end));
for i=1:2
    path(pos_freq-10:pos_freq-9) = num2str(rem(No,8)+i-1,'%02d');
    path_name1(1:2) = num2str(rem(No,8)+i-1,'%02d');
    for j=1:3
        path(end-2:end-1)= num2str(posA+(j-1)*16,'%02d');
        path_name2(8:9)= num2str(posA+(j-1)*16,'%02d');
        file(3:4) = num2str(posB,'%02d');
        file(9:10) = num2str(posA+(j-1)*16,'%02d');
        
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


        % 時間窓設定
        idxmbase=fix((1.5*192000/T-50+0.18*192))+1;
        if idxmbase < 0
            idxmbase = 1;
        end
        time.sppos(i,j)=idxmbase;
        time.sp(i,j)=xtime(idxmbase);
        time.ep(i,j)=xtime(idxmbase+timew);

        % 振幅最大値
        [tymax,tpos]=max(abs(wdm(idxmbase:idxmbase+timew)));
        tpos=idxmbase+tpos;
        tmax.tppos(i,j)=tpos;
        tmax.tpy(i,j)=tymax;
        tmax.tp(i,j)=xtime(tpos);

        % 伝搬時間
        pt.tppos(i,j)=find(abs(wdm(idxmbase:idxmbase+timew))>(tymax*0.1), 1 )+idxmbase;
        pt.tp(i,j)=xtime(pt.tppos(i,j))-0.18;
        
        timesub.sppos(i,j)=pt.tppos(i,j)+fix(0.5*192);
        timesub.sp(i,j)=xtime(timesub.sppos(i,j))-0.18;
        timesub.eppos(i,j)=timesub.sppos(i,j)+fix(0.1*192);
        
        amp=0;
        for k=timesub.sppos(i,j):timesub.eppos(i,j)
            amp=amp+abs(wdm(k));
        end
        amp=amp/(timesub.eppos(i,j)-timesub.sppos(i,j)+1);
        Amp(i,j)=amp;
       
        
        figure(1);subplot(1,2,1)
        hold on
        plot(pt.tp(i,j),Amp(i,j),'o','MarkerEdgeColor','red')
        hold off
        xlim([4.4 4.45]);ylim([0 1.2]);
        xlabel('Time(ms)')
        ylabel('Amplitude')
        grid on
        tp=title({strcat(num2str(rem(No,8)),path_name1(3:end));strcat(num2str(rem(No,8)+1),path_name1(3:end))});set(tp,tx)
    
    end
end


%% from templeture to speed of sound
function T=tem1(t)
    T=331.5+0.6*t;
end


function tp=denpan(path,file,timew)
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


        % 時間窓設定
        idxmbase=fix((1.5*192000/T-50+0.18*192))+1;
        if idxmbase < 0
            idxmbase = 1;
        end

        % 振幅最大値
        [tymax]=max(abs(wdm(idxmbase:idxmbase+timew)));
  
        % 伝搬時間
        tppos=find(abs(wdm(idxmbase:idxmbase+timew))>(tymax*0.1), 1 )+idxmbase;
        tp=xtime(tppos)-0.18;
end