clc
clear
%% ��С���˷� ��һ���Ծ������Ϊ���Ĵ��ۺ����� ������Ϸ���

%������Ϸ���  ��������/������

% С��ʱ��(xi)��λ�ƹ�ϵ(yi)��ϵ
x = [0 1 2 3 4 5 6  7  8  9];
y = [0 2 4 7 8 9 12 14 15 18];

%{
    subplot(m,n,p) ����ǰ�������� m��n��ָ�����ͼ�ֳ� m*n��դ��
    ÿ��դ���� p ����ţ�������ǰ��У����ţ���ŵģ����ԣ��� m = 2��n = 2ʱ��Ź���Ϊ

        1 | 2
        ------
        3 | 4

    ����subplot(2,2,[1 3])����˵������һ����ͼռ�ݵ��� 1�� 3����դ��
    ��subplot(2,2,2)˵����ͼ��ռ�ݵ�2��դ��.
%}
subplot(1,2,1);
plot(x,y,'o');
% ͼ�ε�һЩ����
xlabel('ʱ�䣨�룩');
ylabel('λ�ƣ��ף�');
title('ԭʼ������ɢ��')  
%{
    grid on���Ǵ�����
    grid off���ǹر�����
    ��grid���л�����״̬,�����grid off��״̬��,����grid,�൱��grid on
    �෴,�����grid on״̬������grid �ȼ���grid off
%}
grid on

%{
    polyfit������matlab�����ڽ���������ϵ�һ������������ѧ��������С���˷��������ԭ����
    ������ϣ���֪��ɢ���ϵ����ݼ�������֪�ڵ㼯�ϵĺ���ֵ������һ��������������ͼ��Ϊһ���ߣ�ʹ��ԭ��ɢ���Ͼ����ܽӽ�������ֵ��
    ���÷�����polyfit(x,y,n)���ö���ʽ�����֪��ı���ʽ��
        ����xΪԴ���ݵ��Ӧ�ĺ����꣬��Ϊ������������
            yΪԴ���ݵ��Ӧ�������꣬��Ϊ������������
            nΪ��Ҫ��ϵĽ�����һ��ֱ����ϣ�������������ϣ����ǽ״�Խ��Խ�ã���������������

    ����ʽ��x����ֵy�������������㣺y=polyval(a,x,m)
%}
p = polyfit(x,y,1)
% 0:0.01:9    ��ʼΪ0���յ�Ϊ9������0.01
x1 = 0:0.01:9;
y1 = polyval(p,x1);

x2 = 0:0.01:9;
%{
    MATLAB�еĲ�ֵ����Ϊinterp1������ø�ʽΪ��  yi= interp1(x,y,xi,'method')           
    ����x��yΪ��ֵ�㣬yiΪ�ڱ���ֵ��xi���Ĳ�ֵ�����x,yΪ������ 
    'method'��ʾ���õĲ�ֵ������MATLAB�ṩ�Ĳ�ֵ�����м��֣� 
        'nearest'�����ڽ���ֵ�� 'linear'���Բ�ֵ�� 'spline'����������ֵ�� 'pchip'������ֵ��ȱʡʱ��ʾ���Բ�ֵ
    ע�⣺���еĲ�ֵ������Ҫ��x�ǵ����ģ�����xi���ܹ�����x�ķ�Χ��
%}
y2 = interp1(x,y,x2,'spline');
subplot(1,2,2);
plot(x1,y1,'k',x2,y2,'r')
xlabel('ʱ�䣨�룩');
ylabel('λ�ƣ��ף�');
title('����Ϊ��С���˷���ϣ���ɫΪ��ֵ�����')  
grid on
