function fx=x13(b,x)
%fy=(exp((-b(1)).*x))./(b(2)+b(3).*x);
fx=b(1)*exp(-b(2)*(x-b(3)).^2)+b(4)*exp(-b(5)*(x-b(6)).^2)+b(7)*exp(-b(8)*(x-b(9)).^2)+b(10)*exp(-b(11)*(x-b(12)).^2)+b(13)*exp(-b(14)*(x-b(15)).^2);