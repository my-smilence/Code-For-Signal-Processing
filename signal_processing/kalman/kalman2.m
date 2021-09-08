%% kalman 温度测量
% X(k)=AX(k-1)+W(k)  Z(k)=H*X(k)+V(k);
% X(k)=X(k-1)+K(Z(k-1)-HX(k-1));
% K=P_*（H*p*H'+R）^-1
% P_=F*P*F'+Q
clc
clear

N=120;
CON=25; %房间温度在25摄氏度左右,初始值
Xexpect=CON*ones(1,N);  %初始化范围矩阵 此为期望值E(X)
X=zeros(1,N);   
Xkf=zeros(1,N);     
Z=zeros(1,N);   
P=zeros(1,N);   

X(1)=25.1; %初值 25.1
P(1)=0.01; %P
Z(1)=24.9;
Xkf(1)=Z(1); %

Q=0.01;%W(k)的方差
R=0.25;%V(k)的方差
W=sqrt(Q)*randn(1,N);
V=sqrt(R)*randn(1,N);
F=1;
G=1;
H=1;
I=eye(1);
%%
%%%%%%%%%%%%%%%%%%%%%%%
for k=2:N
    X(k)=F*X(k-1)+G*W(k-1);  
    Z(k)=H*X(k)+V(k);
    X_pre=F*Xkf(k-1);    %状态方程       
    P_pre=F*P(k-1)*F'+Q;  %P更新
    
    Kg=P_pre*inv(H*P_pre*H'+R);  %K增益
    e=Z(k)-H*X_pre;            
    Xkf(k)=X_pre+Kg*e;         
    P(k)=(I-Kg*H)*P_pre;    
end
%%%%%%%%%%%%%%%%
Err_Messure=zeros(1,N);
Err_Kalman=zeros(1,N);
for k=1:N
    Err_Messure(k)=abs(Z(k)-X(k));
    Err_Kalman(k)=abs(Xkf(k)-X(k));
end
t=1:N;
figure('Name','Kalman Filter Simulation','NumberTitle','off');
plot(t,Xexpect,'-b',t,X,'-r',t,Z,'-k',t,Xkf,'-g');
legend('expected','real','measure','kalman extimate');         
xlabel('sample time');
ylabel('temperature');
title('Kalman Filter Simulation');
figure('Name','Error Analysis','NumberTitle','off');
plot(t,Err_Messure,'-b',t,Err_Kalman,'-k');
legend('messure error','kalman error');         
xlabel('sample time');
%%%%%%%%%%%%%%%%%%%%%%%%%

