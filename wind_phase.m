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
pos = str2double(path_name2(8:9));
pos = rem(pos,16);
if pos == 0
    pos = 16;
end
micpos = [12;11;10;9;16;15;14;13];


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

% 時間窓設定
idxmbase=fix(1.5*192000/T+0.18*192);
if idxmbase < 0
idxmbase = 1;
end
time.s=xtime(fix(idxmbase-0.18*192));
time.e=xtime(fix(idxmbase-0.18*192+timew));

wdmpoint = fix(idxmbase+96);


for i=1:2
    path(pos_freq-10:pos_freq-9) = num2str(No+i-1,'%02d');
    path_name1(1:2) = num2str(No+i-1,'%02d');
    %% 2往路
    file(3:4)=num2str(micpos(pos),'%02d');
    for j=1:3
        path(end-2:end-1)= num2str(pos+(j-1)*16,'%02d');
        path_name2(8:9)= num2str(pos+(j-1)*16,'%02d');
        file(9:10) = num2str(pos+(j-1)*16,'%02d');
        
        wdmav=denpan(iftr,d,delay,Fc,path,file,wdmpoint);
        
        figure(1);subplot(1,2,2)
        hold on
        pl = plot(wdmav,'o','MarkerEdgeColor','black'); set(pl,ps)
        set(gca,ax)
        xlim([-1 1]);ylim([-1 1])
        plot([-1 1],[0 0],'k')
        plot([0 0],[-1 1],'k')
        tp=title({strcat(num2str(No),path_name1(3:end));strcat(num2str(No+1),path_name1(3:end))});set(tp,tx)
    end
    %% 3往路
    file(3:4)=num2str(micpos(pos+1),'%02d');
    for j=1:3
        path(end-2:end-1)= num2str(pos+1+(j-1)*16,'%02d');
        path_name2(8:9)= num2str(pos+1+(j-1)*16,'%02d');
        file(9:10) = num2str(pos+1+(j-1)*16,'%02d');

        wdmav=denpan(iftr,d,delay,Fc,path,file,wdmpoint);

        figure(1);subplot(1,2,2)
        plot(wdmav,'o','MarkerEdgeColor','red')
    end

    %% 2復路
    file(3:4)=num2str(pos,'%02d');
    for j=1:3
        path(end-2:end-1)= num2str(micpos(pos)+(j-1)*16,'%02d');
        path_name2(8:9)= num2str(micpos(pos)+(j-1)*16,'%02d');
        file(9:10) = num2str(micpos(pos)+(j-1)*16,'%02d');
        
        wdmav=denpan(iftr,d,delay,Fc,path,file,wdmpoint);
        
        figure(1);subplot(1,2,2)
        plot(wdmav,'s','MarkerEdgeColor','black')
    end
    %% 3復路
    file(3:4)=num2str(pos+1,'%02d');
    for j=1:3
        path(end-2:end-1)= num2str(micpos(pos+1)+(j-1)*16,'%02d');
        path_name2(8:9)= num2str(micpos(pos+1)+(j-1)*16,'%02d');
        file(9:10) = num2str(micpos(pos+1)+(j-1)*16,'%02d');
        
        wdmav=denpan(iftr,d,delay,Fc,path,file,wdmpoint);
        
        figure(1);subplot(1,2,2)
        plot(wdmav,'s','MarkerEdgeColor','red')
    end    
    
    %% 6往路
    file(3:4) = num2str(micpos(8-pos),'%02d');
    for j=1:3
        path(end-2:end-1)= num2str(8-pos+(j-1)*16,'%02d');
        path_name2(8:9)= num2str(8-pos+(j-1)*16,'%02d');
        file(9:10) = num2str(8-pos+(j-1)*16,'%02d');
        wdmav=denpan(iftr,d,delay,Fc,path,file,wdmpoint);                
        figure(1);subplot(1,2,2)
        plot(wdmav,'o','MarkerEdgeColor','blue')
    end
    %% 7往路
    file(3:4) = num2str(micpos(9-pos),'%02d');
    for j=1:3
        path(end-2:end-1)= num2str(9-pos+(j-1)*16,'%02d');
        path_name2(8:9)= num2str(9-pos+(j-1)*16,'%02d');
        file(9:10) = num2str(9-pos+(j-1)*16,'%02d');
        wdmav=denpan(iftr,d,delay,Fc,path,file,wdmpoint);        
        figure(1);subplot(1,2,2)
        plot(wdmav,'o','MarkerEdgeColor','magenta')

    end    
    
    %% 6復路
    file(3:4)=num2str(8-pos,'%02d');
    for j=1:3
        path(end-2:end-1)= num2str(micpos(8-pos)+(j-1)*16,'%02d');
        path_name2(8:9)= num2str(micpos(8-pos)+(j-1)*16,'%02d');
        file(9:10) = num2str(micpos(8-pos)+(j-1)*16,'%02d');
        
        wdmav=denpan(iftr,d,delay,Fc,path,file,wdmpoint);        
        figure(1);subplot(1,2,2)
        plot(wdmav,'s','MarkerEdgeColor','blue')
    end
    %% 7復路
    file(3:4)=num2str(9-pos,'%02d');
    for j=1:3
        path(end-2:end-1)= num2str(micpos(9-pos)+(j-1)*16,'%02d');
        path_name2(8:9)= num2str(micpos(9-pos)+(j-1)*16,'%02d');
        file(9:10) = num2str(micpos(9-pos)+(j-1)*16,'%02d');
        
        wdmav=denpan(iftr,d,delay,Fc,path,file,wdmpoint);        
        figure(1);subplot(1,2,2)
        plot(wdmav,'s','MarkerEdgeColor','magenta')
    end     
end

figure(1);subplot(1,2,2)
dim = [.6 .5 .8 .3 ];
str1=' 2 \rightarrow 11 ';
str2=' 3 \rightarrow 10 ';
str3=' 11 \rightarrow 2 ';
str4=' 10 \rightarrow 3 ';
str5=' 6 \rightarrow 15 ';
str6=' 7 \rightarrow 14 ';
str7=' 15 \rightarrow 6 ';
str8=' 14 \rightarrow 7 ';
str={str1,str2,str3,str4,str5,str6,str7,str8};
annotation('textbox',dim,'String',str,'FitBoxToText','on');

%% デフォルト
path = erase(path,path(pos_wn+4:pos_frame-2));
subpos_wn = strfind(path_name1,'wn');
path_name1 = erase(path_name1,path_name1(subpos_wn+4:end));
for i=1:2
    path(pos_freq-10:pos_freq-9) = num2str(i,'%02d');
    path_name1(1:2) = num2str(i,'%02d');
    file(3:4)=num2str(micpos(pos),'%02d');
    for j=1:3
        path(end-2:end-1)= num2str(pos+(j-1)*16,'%02d');
        path_name2(8:9)= num2str(pos+(j-1)*16,'%02d');
        file(9:10) = num2str(pos+(j-1)*16,'%02d');
                
        wdmav=denpan(iftr,d,delay,Fc,path,file,wdmpoint);
        
        figure(1);subplot(1,2,1)
        hold on
        pl = plot(wdmav,'o','MarkerEdgeColor','black'); set(pl,ps)
        set(gca,ax)
        xlim([-1 1]);ylim([-1 1])
        plot([-1 1],[0 0],'k')
        plot([0 0],[-1 1],'k')
        tp=title({strcat(num2str(rem(No,8)),path_name1(3:end));strcat(num2str(rem(No+1,8)),path_name1(3:end))});set(tp,tx)        
   end
   file(3:4)=num2str(micpos(pos+1),'%02d');
   for j=1:3
        path(end-2:end-1)= num2str(pos+1+(j-1)*16,'%02d');
        path_name2(8:9)= num2str(pos+1+(j-1)*16,'%02d');
        file(9:10) = num2str(pos+1+(j-1)*16,'%02d');
        
        wdmav=denpan(iftr,d,delay,Fc,path,file,wdmpoint);
        
        figure(1);subplot(1,2,1)
        plot(wdmav,'o','MarkerEdgeColor','red')
   end
     
    %% 反対側   
    file(3:4)=num2str(pos,'%02d');
    for j=1:3
        path(end-2:end-1)= num2str(micpos(pos)+(j-1)*16,'%02d');
        path_name2(8:9)= num2str(micpos(pos)+(j-1)*16,'%02d');
        file(9:10) = num2str(micpos(pos)+(j-1)*16,'%02d');
        
        wdmav=denpan(iftr,d,delay,Fc,path,file,wdmpoint);        
        figure(1);subplot(1,2,1)
        plot(wdmav,'s','MarkerEdgeColor','black')
    end
    file(3:4)=num2str(pos+1,'%02d');
    for j=1:3
        path(end-2:end-1)= num2str(micpos(pos+1)+(j-1)*16,'%02d');
        path_name2(8:9)= num2str(micpos(pos+1)+(j-1)*16,'%02d');
        file(9:10) = num2str(micpos(pos+1)+(j-1)*16,'%02d');
        
        wdmav=denpan(iftr,d,delay,Fc,path,file,wdmpoint);        
        figure(1);subplot(1,2,1)
        plot(wdmav,'s','MarkerEdgeColor','red')
    end    
    
    %% 比較対象
    file(3:4) = num2str(micpos(8-pos),'%02d');
    for j=1:3
        path(end-2:end-1)= num2str(8-pos+(j-1)*16,'%02d');
        path_name2(8:9)= num2str(8-pos+(j-1)*16,'%02d');
        file(9:10) = num2str(8-pos+(j-1)*16,'%02d');
        wdmav=denpan(iftr,d,delay,Fc,path,file,wdmpoint);        
        figure(1);subplot(1,2,1)
        plot(wdmav,'o','MarkerEdgeColor','blue')
    end
    
    file(3:4) = num2str(micpos(9-pos),'%02d');
    for j=1:3
        path(end-2:end-1)= num2str(9-pos+(j-1)*16,'%02d');
        path_name2(8:9)= num2str(9-pos+(j-1)*16,'%02d');
        file(9:10) = num2str(9-pos+(j-1)*16,'%02d');
        wdmav=denpan(iftr,d,delay,Fc,path,file,wdmpoint);        
        figure(1);subplot(1,2,1)
        plot(wdmav,'o','MarkerEdgeColor','green')

    end    
    
    %% 反対側
    file(3:4)=num2str(8-pos,'%02d');
    for j=1:3
        path(end-2:end-1)= num2str(micpos(8-pos)+(j-1)*16,'%02d');
        path_name2(8:9)= num2str(micpos(8-pos)+(j-1)*16,'%02d');
        file(9:10) = num2str(micpos(8-pos)+(j-1)*16,'%02d');
        
        wdmav=denpan(iftr,d,delay,Fc,path,file,wdmpoint);        
        figure(1);subplot(1,2,1)
        plot(wdmav,'s','MarkerEdgeColor','blue')
    end
    file(3:4)=num2str(9-pos,'%02d');
    for j=1:3
        path(end-2:end-1)= num2str(micpos(9-pos)+(j-1)*16,'%02d');
        path_name2(8:9)= num2str(micpos(9-pos)+(j-1)*16,'%02d');
        file(9:10) = num2str(micpos(9-pos)+(j-1)*16,'%02d');
        
        wdmav=denpan(iftr,d,delay,Fc,path,file,wdmpoint);        
        figure(1);subplot(1,2,1)
        plot(wdmav,'s','MarkerEdgeColor','green')
    end     
end


%% from templeture to speed of sound
function T=tem1(t)
    T=331.5+0.6*t;
end



function wdmav=denpan(iftr,d,delay,Fc,path,file,wdmpoint)
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
        wdmav=0;
        for k=0:9
        wdmav=wdm(wdmpoint+k)+wdmav;
        end
        wdmav=wdmav/10;
end
