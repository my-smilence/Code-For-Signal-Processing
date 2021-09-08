%%
function [a,b]=Linear_fitting(x,y)
n=size(x,2);
sum_xy=sum(x.*y);
sum_x=sum(x);
sum_y=sum(y);
sum_xx=sum(x.*x);

a=(n*sum_xy-sum_x*sum_y)/(n*sum_xx-sum_x^2);
b=(sum_xy-sum_xx*sum_y/sum_x)/(sum_x-n*sum_xx/sum_x);

% a=(n*sum(x.*y)-sum(x)*sum(y))/(n*sum(x.*x)-sum(x)*sum(x));
% b=(sum(x.*y)-sum(x.*x)*sum(y)/sum(x))/(sum(x)-n*sum(x.*x)/sum(x));
end
