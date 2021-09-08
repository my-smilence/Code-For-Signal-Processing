function Polynomial_fitting_plot(p,t,m)
z=zeros(1,size(t,2));
for i=0:m
    z=z+p(m-i+1)*t.^i;
end
plot(t,z);
end