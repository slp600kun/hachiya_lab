filename = '6_sonar_tr4x4_pos_b41.csv';
T = readtable(filename);
S = table2struct(T);

l1=zeros(1,16);
for i=1:16
    l1(1,i)=leng(S(1).xpos,S(1).ypos,S(i).xpos,S(i).ypos);
end

function len=leng(x1,y1,x2,y2)
    l=(x1-x2).^2+(y1-y2).^2;
    len=sqrt(l);
end



