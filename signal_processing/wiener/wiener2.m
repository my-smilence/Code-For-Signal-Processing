%% ά���˲���
%Ҫ�����ƽ���������
%��ֵ�ź���֪������
%�弤��Ӧ������غ���+����غ�������  Rxx*H=Rxy  H=inv(Rxx)*Rxy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�����ź�
A=1;                                                      %�źŵķ�ֵ
f=1000;                                                 %�źŵ�Ƶ��
fs=10^5;                                                %����Ƶ��
t=(0:999);                                              %������
Mlag=100;                                             %��غ������ȱ���   
x=A*cos(2*pi*f*t/fs);                                %�������Ҳ��ź�
xmean=mean(x);                                    %���Ҳ��źž�ֵ
xvar=var(x,1);                                         %���Ҳ��źŷ���
noise=wgn(1,1000,2);%����1��1000�еľ���ǿ��Ϊ2dbw
xn=x+noise;                                         %�����Ҳ��źż��������Ϊ20dB�ĸ�˹������
plot(t,xn)    
xlabel('x�ᵥλ��t/s','color','b')
ylabel('y�ᵥλ��A/V','color','b')
xnmean = mean(xn)   ;                               %��������źž�ֵ
xnms = mean(xn.^2)    ;                              %��������źž���ֵ
xnvar = var(xn,1)    ;                                   %���������źŷ���
Rxn=xcorr(xn,Mlag,'biased');                   %��������ź�����غ���


%%
figure(2)
subplot(221)
plot((-Mlag:Mlag),Rxn)                             %��������غ���ͼ��
title('�����ź�����غ���ͼ��')
[f,xi]=ksdensity(xn);                                  %��������źŵĸ����ܶȣ�fΪ������xi���ĸ����ܶ�
subplot(222)
plot(xi,f)                                                   %���Ƹ����ܶ�ͼ��
title('�����źŸ����ܶ�ͼ��')
X=fft(xn);                                                  %��������ź����еĿ�����ɢ����Ҷ�任
Px=X.*conj(X)/600;                                   %�����ź�Ƶ��
subplot(223)
semilogy(t,Px)                                          %�����ڰ��������ϵ��Ƶ��ͼ��
title('�����ź��ڰ��������ϵ��Ƶ��ͼ��')
xlabel('x�ᵥλ��w/rad','color','b')
ylabel('y�ᵥλ��w/HZ','color','b')
pxx=periodogram(xn);                               %��������źŵĹ������ܶ�
subplot(224)
semilogy(pxx)                                           %�����ڰ��������ϵ�¹������ܶ�ͼ��
title('�����ź��ڰ��������ϵ�¹������ܶ�ͼ��')
 
xlabel('x�ᵥλ��w/rad','color','b')
ylabel('y�ᵥλ��w/HZ','color','b')
%% 
%ά���˲�
N=100;                                                        %ά���˲�������
Rxnx=xcorr(xn,x,Mlag,'biased');                   %���������ź���ԭʼ�źŵĻ���غ���
rxnx=zeros(N,1);                                       
rxnx(:)=Rxnx(101:101+N-1);
Rxx=zeros(N,N);                                          %���������ź�����ؾ���
Rxx=diag(Rxn(101)*ones(1,N));
for i=2:N
    c=Rxn(101+i)*ones(1,N+1-i);
    Rxx=Rxx+diag(c,i-1)+diag(c,-i+1);
end
Rxx;
h=zeros(N,1);
h=inv(Rxx)*rxnx;                                          %����ά���˲�����h(n)
yn=filter(h,1,xn);                                         %�������ź�ͨ��ά���˲���

%%
figure(5)
plot(yn)                                                      %���ƾ���ά���˲������ź�ͼ��
title('����ά���˲������ź��ź�ͼ��')
xlabel('x�ᵥλ��f/HZ','color','b')
ylabel('y�ᵥλ��A/V','color','b')
ynmean=mean(yn)                                     %���㾭��ά���˲������źž�ֵ
ynms=mean(yn.^2)                                     %���㾭��ά���˲������źž���ֵ
ynvar=var(yn,1)                                         %���㾭��ά���˲������źŷ���
Ryn=xcorr(yn,Mlag,'biased');                     %���㾭��ά���˲������ź�����غ���
figure(6)
subplot(221)
plot((-Mlag:Mlag),Ryn)                               %��������غ���ͼ��
title('����ά���˲������ź�����غ���ͼ��')
[f,yi]=ksdensity(yn);                                    %���㾭��ά���˲������źŵĸ����ܶȣ�fΪ������xi���ĸ����ܶ�
subplot(222)
plot(yi,f)                                                     %���Ƹ����ܶ�ͼ��
title('����ά���˲������źŸ����ܶ�ͼ��')
Y=fft(yn);                                                   %���㾭��ά���˲������ź����еĿ�����ɢ����Ҷ�任
Py=Y.*conj(Y)/600;                                    %�����ź�Ƶ��
subplot(223)
semilogy(t,Py)                                           %�����ڰ��������ϵ��Ƶ��ͼ��
title('����ά���˲������ź��ڰ��������ϵ��Ƶ��ͼ��')
xlabel('x�ᵥλ��w/rad','color','b')
ylabel('y�ᵥλ��w/HZ','color','b')
pyn=periodogram(yn);                               %���㾭��ά���˲������źŵĹ������ܶ�
subplot(224)
semilogy(pyn)                                            %�����ڰ��������ϵ�¹������ܶ�ͼ��
title('����ά���˲������ź��ڰ��������ϵ�¹������ܶ�ͼ��')
xlabel('x�ᵥλ��w/rad','color','b')
ylabel('y�ᵥλ��w/HZ','color','b')
subplot(4,1,1),plot(noise); title('�����ź�')
subplot(4,1,2),plot(x); title('�����ź�')
subplot(4,1,3),plot(xn); title('�����ź�')
subplot(4,1,4),plot(yn); title('ά���ź�')


