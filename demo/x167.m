function fx=x167(b,x)
%fx=b(1)+b(2)*x+b(3)*x.^2./1+b(4)*x;
fx=b(1).*exp(-(x./b(2)))+b(3).*exp(-(x./b(4)))+b(5);