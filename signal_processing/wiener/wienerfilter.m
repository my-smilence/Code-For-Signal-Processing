%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%维纳滤波器函数设计
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y =wienerfilter(x,Rxx,Rxd,N) 
%进行维纳滤波 
%x是输入信号，Rxx是输入信号的自相关向量 
%Rxd是输入信号和理想信号的的互相关向量，N是维纳滤波器的长度 
%输出y是输入信号通过维纳滤波器进行维纳滤波后的输出 
h=yulewalker(Rxx,Rxd,N);								%求解维纳滤波器系数 
t=conv(x,h);											%进行滤波 
Lh=length(h);											%得到滤波器的长度 
Lx=length(x);											%得到输入信号的长度 
y=t(double(uint16(Lh/2)):Lx+double(uint16(Lh/2))-1);%输出序列y的长度和输入序列x的长度相同
%以下是维纳滤波器系数的求解 
function h=yulewalker(A,B,M)    
%求解Yule-Walker方程 
%A是接收信号的自相关向量为 Rxx(0),Rxx(1),......,Rxx(M-1) 
%B是接收信号和没有噪声干扰信号的互相关向量为 Rxd(0),Rxd(1),......,Rxd(M-1) 
%M是滤波器的长度 
%h保存滤波器的系数 
T1=zeros(1,M);%T1存放中间方程的解向量 
T2=zeros(1,M);%T2存放中间方程的解向量 
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
%维纳滤波器案例测试
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear all; close all;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load handel						%加载语音信号
d=y; d=d*8;						%增强语音信号强度
d=d';
[m,n]=size(d);
T = 1/Fs; % 采样时间
t = (1:n)*T;% 时间
subplot(3,2,1);
plot(t,d);				
title('原始语音信号');
xlabel('时间/t');
ylabel('幅值/dB');
fq=fft(d,8192);						%进行傅立叶变换得到语音信号频频
subplot(3,2,2);
f=Fs*(0:4095)/8192;
plot(f,abs(fq(1:4096)));				%画出频谱图
title('原始语音信号的频域图形');
xlabel('频率 f');
ylabel('FFT');
x_noise=randn(1,n);				%（0，1）分布的高斯白噪声
x=d+x_noise;						%加入噪声后的语音信号
subplot(3,2,3);
plot(t,x);				
title('加入噪声后');
xlabel('时间/t');
ylabel('幅值/dB');
fq=fft(x,8192);						%对加入噪声后的信号进行傅立叶变换，看其频谱变化
subplot(3,2,4);
plot(f,abs(fq(1:4096)));				%画出加入噪声后信号的频谱图
title('加入噪声后语音信号的频域图形');
xlabel('频率 f');
ylabel('FFT');
%维纳滤波
yyhxcorr=xcorr(x(1:4096));			%求取信号的信号的自相关函数
size(yyhxcorr); 
A=yyhxcorr(4096:4595);
yyhdcorr=xcorr(d(1:4096),x(1:4096));			%求取信号和理想信号的互相关函数
size(yyhdcorr);
B=yyhdcorr(4096:4595);
M=500;
yyhresult=wienerfilter(x,A,B,M);				%进行维纳滤波
yyhresult=yyhresult(300:8192+299);
subplot(3,2,5);
t = (1:8192)*T;% 时间
plot(t,yyhresult);				%画出频谱图
title('进行维纳滤波');
xlabel('时间/t');
ylabel('幅值/dB');
fq=fft(yyhresult);							%对维纳滤波的结果进行傅立叶变换，看其频谱变化
subplot(3,2,6); 
f=Fs*(0:4095)/8192;
plot(f,abs(fq(1:4096)));						%画出维纳滤波后信号的频谱图
title('经过维纳滤波后语音信号的频域图形');
xlabel('频率 f');
ylabel('FFT');
