%% 2-----���ݷ���
detd5=2048*1;
i5=0;
dtf5=i5;    %���ݴ��ڣ����ȣ��Ļ���λ��
tf5=(dtf5+1:1:dtf5+detd5);  %���������ݳ��ȣ�

numPart=20;
%����APD���յ��źŵ�һ��
APDpart1=zeros(detd5,numPart*(numData+numDataI));
for i=1:1:(numData+numDataI)
    for j=1:1:numPart
        APDpart1(:,(i-1)*numPart+j)=APD( (detd5*(j-1)+1):detd5*j,i);
    end
end
%���������źŵ�һ��
tzpart1=zeros(detd5,numPart*(numData+numDataI));
for i=1:1:(numData+numDataI)
    for j=1:1:numPart
        tzpart1(:,(i-1)*numPart+j)=tz( (detd5*(j-1)+1):detd5*j,i);
    end
end
%���������źŵ�һ��
jidaipart1=zeros(detd5,numPart*(numData+numDataI));
for i=1:1:(numData+numDataI)
    for j=1:1:numPart
        jidaipart1(:,(i-1)*numPart+j)=jidai( (detd5*(j-1)+1):detd5*j,i);
    end
end
disp('2nd: End of group.');

