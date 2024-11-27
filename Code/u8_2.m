%% 2-----数据分组
detd5=2048*1;
i5=0;
dtf5=i5;    %数据窗口（长度）的滑动位置
tf5=(dtf5+1:1:dtf5+detd5);  %采样的数据长度，

numPart=20;
%读出APD接收到信号的一份
APDpart1=zeros(detd5,numPart*(numData+numDataI));
for i=1:1:(numData+numDataI)
    for j=1:1:numPart
        APDpart1(:,(i-1)*numPart+j)=APD( (detd5*(j-1)+1):detd5*j,i);
    end
end
%读出调制信号的一份
tzpart1=zeros(detd5,numPart*(numData+numDataI));
for i=1:1:(numData+numDataI)
    for j=1:1:numPart
        tzpart1(:,(i-1)*numPart+j)=tz( (detd5*(j-1)+1):detd5*j,i);
    end
end
%读出基带信号的一份
jidaipart1=zeros(detd5,numPart*(numData+numDataI));
for i=1:1:(numData+numDataI)
    for j=1:1:numPart
        jidaipart1(:,(i-1)*numPart+j)=jidai( (detd5*(j-1)+1):detd5*j,i);
    end
end
disp('2nd: End of group.');

