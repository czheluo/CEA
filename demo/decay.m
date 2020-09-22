clear,clc
%warning('off',msgID);
x=xlsread('decay','0d','a2:a36');
y=xlsread('decay','0d','s2:s36');
n=length(y);%%
%fx=@(b,x)(b(1).*(x.^2+b(2).*x))./(x.^2+b(3).*x+b(4));
fx=@(b,x)b(1).*exp(-(x./b(2)))+b(3).*exp(-(x./b(4)))+b(5);
%b=[ 13.703       1.2713      0.04604     0.094421];
b=[69.65856712  0.267654987  349.0694531  0.01792381726  22.41263325];
%fx=@(b,t)b(1)*exp(-t/b(2))+b(3)*exp(-t/b(4))+b(5);
%b=[257.2038649  1.908018941e-05  33.5724231  0.0002621509447  14.3301809];
%b=[257.2      0.01908       33.572      0.26215        14.33];
for l=1:20
    b=lsqcurvefit(fx,b,x,y)
    b=nlinfit(x,y,fx,b)
end
figure(1),clf
plot(x,y,'o','markersize',10,'markerfacecolor','k')
hold on
x1=linspace(min(x)-range(x)/150,max(x)+range(x)/50,150);
y1=fx(b,x1);
plot(x1,y1,'r-','linewidth',2.5)
%xlabel('X'),ylabel('Y')
axis tight
legend('data','fitiing','location','best')
SSy=var(y)*(n-1)
y2=fx(b,x);
SSe=(y-y2)'*(y-y2)
r2=(SSy-SSe)/SSy
MSe=SSe/(n-5)
[R,P]=corr(x,y)