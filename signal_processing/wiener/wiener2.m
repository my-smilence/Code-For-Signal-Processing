%% 维纳滤波：
%要求处理宽平稳随机过程
%真值信号已知，先验
%冲激响应由自相关函数+互相关函数决定  Rxx*H=Rxy  H=inv(Rxx)*Rxy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%输入信号
A=1;                                                      %信号的幅值
f=1000;                                                 %信号的频率
fs=10^5;                                                %采样频率
t=(0:999);                                              %采样点
Mlag=100;                                             %相关函数长度变量   
x=A*cos(2*pi*f*t/fs);                                %输入正弦波信号
xmean=mean(x);                                    %正弦波信号均值
xvar=var(x,1);                                         %正弦波信号方差
noise=wgn(1,1000,2);%产生1行1000列的矩阵，强度为2dbw
xn=x+noise;                                         %给正弦波信号加入信噪比为20dB的高斯白噪声
plot(t,xn)    
xlabel('x轴单位：t/s','color','b')
ylabel('y轴单位：A/V','color','b')
xnmean = mean(xn)   ;                               %计算加噪信号均值
xnms = mean(xn.^2)    ;                              %计算加噪信号均方值
xnvar = var(xn,1)    ;                                   %计算输入信号方差
Rxn=xcorr(xn,Mlag,'biased');                   %计算加噪信号自相关函数


%%
figure(2)
subplot(221)
plot((-Mlag:Mlag),Rxn)                             %绘制自相关函数图像
title('加噪信号自相关函数图像')
[f,xi]=ksdensity(xn);                                  %计算加噪信号的概率密度，f为样本点xi处的概率密度
subplot(222)
plot(xi,f)                                                   %绘制概率密度图像
title('加噪信号概率密度图像')
X=fft(xn);                                                  %计算加噪信号序列的快速离散傅里叶变换
Px=X.*conj(X)/600;                                   %计算信号频谱
subplot(223)
semilogy(t,Px)                                          %绘制在半对数坐标系下频谱图像
title('输入信号在半对数坐标系下频谱图像')
xlabel('x轴单位：w/rad','color','b')
ylabel('y轴单位：w/HZ','color','b')
pxx=periodogram(xn);                               %计算加噪信号的功率谱密度
subplot(224)
semilogy(pxx)                                           %绘制在半对数坐标系下功率谱密度图像
title('加噪信号在半对数坐标系下功率谱密度图像')
 
xlabel('x轴单位：w/rad','color','b')
ylabel('y轴单位：w/HZ','color','b')
%% 
%维纳滤波
N=100;                                                        %维纳滤波器长度
Rxnx=xcorr(xn,x,Mlag,'biased');                   %产生加噪信号与原始信号的互相关函数
rxnx=zeros(N,1);                                       
rxnx(:)=Rxnx(101:101+N-1);
Rxx=zeros(N,N);                                          %产生加噪信号自相关矩阵
Rxx=diag(Rxn(101)*ones(1,N));
for i=2:N
    c=Rxn(101+i)*ones(1,N+1-i);
    Rxx=Rxx+diag(c,i-1)+diag(c,-i+1);
end
Rxx;
h=zeros(N,1);
h=inv(Rxx)*rxnx;                                          %计算维纳滤波器的h(n)
yn=filter(h,1,xn);                                         %将加噪信号通过维纳滤波器

%%
figure(5)
plot(yn)                                                      %绘制经过维纳滤波器后信号图像
title('经过维纳滤波器后信号信号图像')
xlabel('x轴单位：f/HZ','color','b')
ylabel('y轴单位：A/V','color','b')
ynmean=mean(yn)                                     %计算经过维纳滤波器后信号均值
ynms=mean(yn.^2)                                     %计算经过维纳滤波器后信号均方值
ynvar=var(yn,1)                                         %计算经过维纳滤波器后信号方差
Ryn=xcorr(yn,Mlag,'biased');                     %计算经过维纳滤波器后信号自相关函数
figure(6)
subplot(221)
plot((-Mlag:Mlag),Ryn)                               %绘制自相关函数图像
title('经过维纳滤波器后信号自相关函数图像')
[f,yi]=ksdensity(yn);                                    %计算经过维纳滤波器后信号的概率密度，f为样本点xi处的概率密度
subplot(222)
plot(yi,f)                                                     %绘制概率密度图像
title('经过维纳滤波器后信号概率密度图像')
Y=fft(yn);                                                   %计算经过维纳滤波器后信号序列的快速离散傅里叶变换
Py=Y.*conj(Y)/600;                                    %计算信号频谱
subplot(223)
semilogy(t,Py)                                           %绘制在半对数坐标系下频谱图像
title('经过维纳滤波器后信号在半对数坐标系下频谱图像')
xlabel('x轴单位：w/rad','color','b')
ylabel('y轴单位：w/HZ','color','b')
pyn=periodogram(yn);                               %计算经过维纳滤波器后信号的功率谱密度
subplot(224)
semilogy(pyn)                                            %绘制在半对数坐标系下功率谱密度图像
title('经过维纳滤波器后信号在半对数坐标系下功率谱密度图像')
xlabel('x轴单位：w/rad','color','b')
ylabel('y轴单位：w/HZ','color','b')
subplot(4,1,1),plot(noise); title('噪声信号')
subplot(4,1,2),plot(x); title('正弦信号')
subplot(4,1,3),plot(xn); title('加噪信号')
subplot(4,1,4),plot(yn); title('维纳信号')


