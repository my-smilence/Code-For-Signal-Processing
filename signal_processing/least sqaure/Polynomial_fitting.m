function P=Polynomial_fitting(x,y,m)  %x,yΪ���г�����ȵ�����������mΪ��϶���ʽ����
a=zeros(2*m+1,1);
for i=0:2*m
    a(i+1,1)=sum(x.^i);
end
b=zeros(m+1,1);
for i=0:m
    b(i+1,1)=sum((x.^i).*y);
end
A=zeros(m+1,m+1);
for i=1:m+1
    for j=1:m+1
        A(i,j)=a(i+j-1,1);
    end
end
p=A\b;
P=fliplr(p'); %���ҷ�ת����
%https://www.cnblogs.com/kailugaji/p/6932482.html
end

function Polynomial_fitting_plot(p,t,m)
z=zeros(1,size(t,2));
for i=0:m
    z=z+p(m-i+1)*t.^i;
end
plot(t,z);
end
