%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ά���˲����������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y =wienerfilter(x,Rxx,Rxd,N) 
%����ά���˲� 
%x�������źţ�Rxx�������źŵ���������� 
%Rxd�������źź������źŵĵĻ����������N��ά���˲����ĳ��� 
%���y�������ź�ͨ��ά���˲�������ά���˲������� 
h=yulewalker(Rxx,Rxd,N);								%���ά���˲���ϵ�� 
t=conv(x,h);											%�����˲� 
Lh=length(h);											%�õ��˲����ĳ��� 
Lx=length(x);											%�õ������źŵĳ��� 
y=t(double(uint16(Lh/2)):Lx+double(uint16(Lh/2))-1);%�������y�ĳ��Ⱥ���������x�ĳ�����ͬ
%������ά���˲���ϵ������� 
function h=yulewalker(A,B,M)    
%���Yule-Walker���� 
%A�ǽ����źŵ����������Ϊ Rxx(0),Rxx(1),......,Rxx(M-1) 
%B�ǽ����źź�û�����������źŵĻ��������Ϊ Rxd(0),Rxd(1),......,Rxd(M-1) 
%M���˲����ĳ��� 
%h�����˲�����ϵ�� 
T1=zeros(1,M);%T1����м䷽�̵Ľ����� 
T2=zeros(1,M);%T2����м䷽�̵Ľ����� 
T1(1)=B(1)/A(1); 
T2(1)=A(2)/A(1); 
X=zeros(1,M); 
for i=2:M-1 
temp1=0; 
temp2=0; 
    for j=1:i-1 
        temp1=temp1+A(i-j+1)*T1(j); 
        temp2=temp2+A(i-j+1)*T2(j); 
    end 
    X(i)=(B(i)-temp1)/(A(1)-temp2); 
    for j=1:i-1 
        X(j)=T1(j)-X(i)*T2(j); 
    end 
    for j=1:i 
        T1(j)=X(j); 
    end 
temp1=0; 
temp2=0; 
    for j=1:i-1 
        temp1=temp1+A(j+1)*T2(j); 
        temp2=temp2+A(j+1)*T2(i-j); 
    end 
    X(1)=(A(i+1)-temp1)/(A(1)-temp2); 
    for j=2:i 
        X(j)=T2(j-1)-X(1)*T2(i-j+1); 
    end 
    for j=1:i 
        T2(j)=X(j); 
    end 
end 
temp1=0; 
temp2=0; 
for j=1:M-1 
	temp1=temp1+A(M-j+1)*T1(j); 
	temp2=temp2+A(M-j+1)*T2(j); 
end 
X(M)=(B(M)-temp1)/(A(1)-temp2); 
for j=1:M-1 
	X(j)=T1(j)-X(M)*T2(j); 
end 
h=X;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ά���˲�����������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear all; close all;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load handel						%���������ź�
d=y; d=d*8;						%��ǿ�����ź�ǿ��
d=d';
[m,n]=size(d);
T = 1/Fs; % ����ʱ��
t = (1:n)*T;% ʱ��
subplot(3,2,1);
plot(t,d);				
title('ԭʼ�����ź�');
xlabel('ʱ��/t');
ylabel('��ֵ/dB');
fq=fft(d,8192);						%���и���Ҷ�任�õ������ź�ƵƵ
subplot(3,2,2);
f=Fs*(0:4095)/8192;
plot(f,abs(fq(1:4096)));				%����Ƶ��ͼ
title('ԭʼ�����źŵ�Ƶ��ͼ��');
xlabel('Ƶ�� f');
ylabel('FFT');
x_noise=randn(1,n);				%��0��1���ֲ��ĸ�˹������
x=d+x_noise;						%����������������ź�
subplot(3,2,3);
plot(t,x);				
title('����������');
xlabel('ʱ��/t');
ylabel('��ֵ/dB');
fq=fft(x,8192);						%�Լ�����������źŽ��и���Ҷ�任������Ƶ�ױ仯
subplot(3,2,4);
plot(f,abs(fq(1:4096)));				%���������������źŵ�Ƶ��ͼ
title('���������������źŵ�Ƶ��ͼ��');
xlabel('Ƶ�� f');
ylabel('FFT');
%ά���˲�
yyhxcorr=xcorr(x(1:4096));			%��ȡ�źŵ��źŵ�����غ���
size(yyhxcorr); 
A=yyhxcorr(4096:4595);
yyhdcorr=xcorr(d(1:4096),x(1:4096));			%��ȡ�źź������źŵĻ���غ���
size(yyhdcorr);
B=yyhdcorr(4096:4595);
M=500;
yyhresult=wienerfilter(x,A,B,M);				%����ά���˲�
yyhresult=yyhresult(300:8192+299);
subplot(3,2,5);
t = (1:8192)*T;% ʱ��
plot(t,yyhresult);				%����Ƶ��ͼ
title('����ά���˲�');
xlabel('ʱ��/t');
ylabel('��ֵ/dB');
fq=fft(yyhresult);							%��ά���˲��Ľ�����и���Ҷ�任������Ƶ�ױ仯
subplot(3,2,6); 
f=Fs*(0:4095)/8192;
plot(f,abs(fq(1:4096)));						%����ά���˲����źŵ�Ƶ��ͼ
title('����ά���˲��������źŵ�Ƶ��ͼ��');
xlabel('Ƶ�� f');
ylabel('FFT');
