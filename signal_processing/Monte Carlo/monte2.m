clear all
%% MCM intrgral for function f=x^2 in (-1,1)
%generating random number in (0,1)
N=100000;
P=rand(N,2);
%extending x and y range linearly
x=2*P(:,1)-1;
y=P(:,2);
%MCM by using 'find' function
R=find(y<=x.^2&y>0);
Num=length(R);
Int=Num/N*2
% plot
figure();
plot(x(R),y(R),'r.')
%annotation
dim = [.2 .5 .3 .3];
result=num2str(Int);
str=['Monte Carlo Method Integration \int_{-1}^{1} x^2 dx = ',result];
annotation('textbox',dim,'String',str,'FitBoxToText','on');