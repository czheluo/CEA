% curve and surface fitting with C-E algorithm
clear all;clc,format short g
global n x x1 x2 x3 x4 y k b b0 bo bi bm db c q1 qm qv t1 s1 n1 nc dim v1 v2 v
!copy x167.m nfx.m;
%X=xlsread('qs','a1:a580');
x=xlsread('instance','1','a2:a15');
y=xlsread('instance','1','b2:b15');
X=[x y];
k=4;
[n,dim]=size(X);
y=X(:,dim);dim=dim-1;
N=ceil(480+25*k^2);
b0=rand(1,k);
b0=[ 691.6188819  50.43652062  0.5040764565 1];
bi=.017*b0;bm=b0;bo=b0;b=b0;
SSy=var(y)*(n-1);
if dim==1
    x=X(:,1);x1a=min(x);x1b=max(x);r1=range(x);
    x11=linspace(x1a-r1/350,x1b+r1/50,350);
elseif dim==2
    x1=X(:,1);x2=X(:,2);X=[];
    x1a=min(x1);x1b=max(x1);r1=range(x1);
    x2a=min(x2);x2b=max(x2);r2=range(x2);
    [x11 x22]=meshgrid(x1a-r1/100:r1/165:x1b,x2a:r2/165:x2b);
elseif dim==3
    x1=X(:,1);x2=X(:,2);x3=X(:,3);
    x1a=min(x1);x1b=max(x1);r1=range(x1);
    x2a=min(x2);x2b=max(x2);r2=range(x2);
    x3a=min(x3);x3b=max(x3);r3=range(x3);
    m1=mean(x1);m2=mean(x2);m3=mean(x3);
elseif dim==4
    x1=X(:,1);x2=X(:,2);x3=X(:,3);x4=X(:,4);
elseif dim==5
    x1=X(:,1);x2=X(:,2);x3=X(:,3);x4=X(:,4);x5=X(:,5);
end
ya=min(y);yb=max(y);ry=yb-ya;ya=ya-ry/30;yb=yb+ry/30;
str=num2str([1:n]');
figure(1),clf
hold off
if dim==1
    y1=nfx(bm,x11);
    plot(x,y,'ko',x11,y1,'b-','markerfacecolor','k','linewidth',3)
    axis tight
    pause(.0001)
elseif dim==2
    plot3(x1,x2,y,'o')
    stem3(x1,x2,y,'filled')
    text(x1,x2,y+.04*ry,str,'fontsize',12)
    pause(.0001)
    hold on
    y1=nfx(bm,x11,x22);
    [np,nq]=size(y1);
    for i=1:np
        for j=1:nq
            if ~isreal(y1(i,j));
                y1(i,j)=nan;
            end
        end
    end
    surf(x11,x22,y1)
    axis([x1a x1b x2a x2b ya yb])
    alpha(.85)
    shading interp
    axis tight
    pause(.0001)
elseif dim==3
    figure(1),clf
    plot3(x1,x2,y,'o')
    stem3(x1,x2,y,'filled')
    text(x1,x2,y+.04*ry,str,'fontsize',10)
    pause(.0001)
    hold on
    [x11,x22]=meshgrid(x1a:r1/50:x1b,x2a:r2/50:x2b);
    y1=nfx(bm,x11,x22,m3);y1=real(y1);
    surf(x11,x22,y1)
    axis([x1a x1b x2a x2b ya yb])
    alpha(.85)
    shading interp
    axis tight
    pause(.0001)
    xlabel('X1'),ylabel('X2'),zlabel('Y')
    figure(2),clf
    [x11,x33]=meshgrid(x1a:r1/50:x1b,x3a:r3/50:x3b);
    plot3(x1,x3,y,'o')
    stem3(x1,x3,y,'filled')
    text(x1,x3,y+.04*ry,str,'fontsize',10)
    pause(.0001)
    hold on
    y2=nfx(bm,x11,m2,x33);
    surf(x11,x33,y2)
    axis([x1a x1b x3a x3b ya yb])
    alpha(.85)
    shading interp
    axis tight
    pause(.0001)
    xlabel('X1'),ylabel('X3'),zlabel('Y')
    figure(3),clf
    plot3(x2,x3,y,'o')
    stem3(x2,x3,y,'filled')
    text(x2,x3,y+.04*ry,str,'fontsize',10)
    pause(.0001)
    hold on
    [x22,x33]=meshgrid(x2a:r2/50:x2b,x3a:r3/50:x3b);
    y3=nfx(bm,m1,x22,x33);
    surf(x22,x33,y3)
    axis([x2a x2b x3a x3b ya yb])
    alpha(.85)
    shading interp
    axis tight
    xlabel('X2'),ylabel('X3'),zlabel('Y')
    pause(.0001)
end
%stop
bmr=zeros(10,k);qmr(1:10)=1e280;bmn=bm;B=zeros(N,k);
t1=zeros(1,k);s1=zeros(1,k);n1=0;nc=0;c=.01*SSy/n+1;
q1=1e280;qm(1:9)=q1;b=bm;sub5,bo=bm;v1=9;qm(v1)=q1;
qmn(1:10)=1e280;v0=10;qmn(v0)=qm(v1);rep=0;v3=1;
while (qmn(v0-3)-qmn(v0))/qmn(v0)>=1e-12 && rep<=3
    for v=1:8+1.4*mod(v1,4)
        for l=1:N
            b=random('unif',b0-1.5*bi,b0+1.5*bi);
            sub5
        end
        b0=bo;bi=bi/2.05;pause(.000001)
    end
    bi=(1.8+4.5e-3*mod(v1,9)^3)*bi;
    nc2=nc;qm3=qm(v1);qm2=qm3+1;v=1;b0=bm;bo=bm;
    while v<=4.3+.8*mod(v1,4) || qm3<qm2
        for l=1:N
            b=random('unif',b0-1.6*bi,b0+1.6*bi);
            sub6
        end
        b0=bo;bi=3.7*bi;v=v+1;qm2=qm3;qm3=qm(v1);
        pause(.000001)
    end
    nc2=nc-nc2;nc1=nc;b0=bm;bo=b0;
    while nc2<=.023*N && c<=1e3*SSy
        c=c*2.3;
        for l=1:N
            b=random('unif',b0-1.5*bi,b0+1.5*bi);
            sub6
        end
        b0=bo;nc2=nc2+nc-nc1;nc1=nc;bi=bi*1.25;
    end
    if isnan(qm(v1))==1 || isreal(bm)==0
        bm=bmn;
        bi=.17*abs(bm);b0=bm;bo=b0;b=b0;
        q1=1e280;t1=zeros(1,k);s1=zeros(1,k);n1=1;nc=1;
        qm(v1)=qmn(v0);
    end
    hold off
    if dim==1
        y11=nfx(bm,x11);
        plot(x,y,'ko',x11,y11,'b-','markerfacecolor','k','linewidth',3)
        axis tight
        pause(.001)
    elseif dim==2
        plot3(x1,x2,y,'o','markerfacecolor','b');
        if n<250
            stem3(x1,x2,y,'filled')
        end
        pause(.01)
        if n<200
            text(x1,x2,y+.04*ry,str,'fontsize',12)
        end
        hold on
        y11=nfx(bm,x11,x22);
        for i=1:np
            for j=1:nq
                if ~isreal(y11(i,j));
                    y11(i,j)=nan;
                end
            end
        end
        surf(x11,x22,y11)
        axis tight
        alpha(.85)
        shading interp
        %y1=nfx(bm,x1,x2);
        %plot3(x1,x2,y1,'r-','linewidth',3)
        pause(.00001)
    elseif dim==3
        figure(1),clf
        plot3(x1,x2,y,'o')
        stem3(x1,x2,y,'filled')
        if n<205
            text(x1,x2,y+.04*ry,str,'fontsize',10)
        end
        pause(.0001)
        hold on
        [x11,x22]=meshgrid(x1a:r1/50:x1b,x2a:r2/50:x2b);
        y1=nfx(bm,x11,x22,m3);
        surf(x11,x22,y1)
        axis([x1a x1b x2a x2b ya yb])
        alpha(.7)
        shading interp
        axis tight
        pause(.0001)
        xlabel('X1'),ylabel('X2'),zlabel('Y')
        figure(2),clf
        [x11,x33]=meshgrid(x1a:r1/50:x1b,x3a:r3/50:x3b);
        plot3(x1,x3,y,'o')
        stem3(x1,x3,y,'filled')
        if n<205
            text(x1,x3,y+.04*ry,str,'fontsize',10)
        end
        pause(.0001)
        hold on
        y2=nfx(bm,x11,m2,x33);
        surf(x11,x33,y2)
        axis([x1a x1b x3a x3b ya yb])
        alpha(.7)
        shading interp
        axis tight
        pause(.0001)
        xlabel('X1'),ylabel('X3'),zlabel('Y')
        figure(3),clf
        plot3(x2,x3,y,'o')
        stem3(x2,x3,y,'filled')
        if n<205
            text(x2,x3,y+.04*ry,str,'fontsize',10)
        end
        pause(.0001)
        hold on
        [x22,x33]=meshgrid(x2a:r2/50:x2b,x3a:r3/50:x3b);
        y3=nfx(bm,m1,x22,x33);
        surf(x22,x33,y3)
        axis([x2a x2b x3a x3b ya yb])
        alpha(.7)
        shading interp
        axis tight
        xlabel('X2'),ylabel('X3'),zlabel('Y')
        pause(.0001)
    end
    v1=v1+1;qm(v1)=qm(v1-1);
    if qm(v1)<=qmn(v0)
        qmn(v0)=qm(v1);
        bmn=bm;
    end
    %disp([qm(v1),bm])
    disp([qm(v1),c,nc,v-1])
    pause(.001)
    if mod(v1,4)==1
        ct=exp(-(nc-.5*N)/N);c=c*ct;
        bi=(.85+.15*exp(1.5-1.5*ct))*bi;
        b=bm;sub5,b0=bm;
    elseif mod(v1,4)==2
        ct=exp(-1.2*(nc-.95*N)/N);c=c*ct;
        if v1>=12 && n1>30
            b=bm;db=1.75e-5*(bi+sqrt((s1-t1.*t1/n1)/n1));
            qv(1:10)=qm(v1);v2=9;
            mqt2
        end
        sb=sqrt((s1-t1.*t1/(n1+(n1==0)))/(n1+(n1==0)));
        bi=.35*bi+1.55*sb;b=.85*bm+.15*t1/(n1+(n1==0));
        sub5,b0=bm;
        s1=zeros(1,k);t1=zeros(1,k);n1=0;nc=0;
    elseif mod(v1,4)==3
        ct=exp(-(nc-.6*N)/N);c=c*ct;
        bi=(.85+.15*exp(1.5-1.5*ct))*bi;
        b=bm;sub5,b0=bm;
    else
        ct=exp(-1.2*(nc-N)/N);c=c*ct;
        if v1>=13 && n1>40
            b=bm;db=1.75e-5*(bi+sqrt((s1-t1.*t1/n1)/n1));
            qv(1:10)=qm(v1);v2=9;
            mqt1
        end
        sb=sqrt((s1-t1.*t1/(n1+(n1==0)))/(n1+(n1==0)));
        bi=.013*bi+(1.8+log(2*qm(v1)/SSy+.999))*sb;% this is the next step length
        b0=.95*bm+.05*t1/(n1+(n1==0));bo=bm;
        s1=zeros(1,k);t1=zeros(1,k);n1=0;nc=0;
    end
    if (qm(v1-5)-qm(v1))/qm(v1)<k*k*1e-12
        if str2double(str2mat(vpa(qm(v1),8)))<str2double(str2mat(vpa(qmn(v0-1),8)))
            rep=0;
        elseif str2double(str2mat(vpa(qm(v1),8)))==str2double(str2mat(vpa(qmn(v0-1),8)))
            rep=rep+1;
        end
        qr=10;
        for i=1:10
            if qm(v1)<=qmr(i)
                qr=i;break
            end
        end
        for i=qr:9
            qmr(i+1)=qmr(i);
            bmr(i+1,:)=bmr(i,:);
        end
        qmr(qr)=qm(v1);bmr(qr,:)=bm;
        b0=.01*rand(1,k)+.35*mod(v0,5)*mean(bmr(1:3,:));
        c=9*c+SSy/n;q1=1e250;qm(1:10)=q1;v1=10;
        for j=1:k
            while abs(b0(j))>1e7
                b0(j)=b0(j)/50;
            end
        end
        bi=abs(b0)*1.5;
        qm(v1)=q1/1.5;bm=b0;bo=b0;
        v0=v0+1;qmn(v0)=qmn(v0-1);
        nc=0;n1=0;t1=zeros(1,k);s1=zeros(1,k);v3=v3+10;
    end
end
SSy
b=sprintf('  %0.10g',bmn),RSS=sprintf(' %0.12g',qmn(v0)),
MSe=qmn(v0)/(n-k-1),R2=(SSy-qmn(v0))/SSy,
b=bmn;
if dim==3
    figure(1),clf
    plot3(x1,x2,y,'o')
    stem3(x1,x2,y,'filled')
    if n<250
        text(x1,x2,y+.04*ry,str,'fontsize',12)
    end
    pause(.0001)
    hold on
    [x11,x22]=meshgrid(x1a:r1/75:x1b,x2a:r2/75:x2b);
    y1=nfx(b,x11,x22,m3);
    surf(x11,x22,y1)
    axis([x1a x1b x2a x2b ya yb])
    alpha(.85)
    shading interp
    axis tight
    pause(.0001)
    xlabel('X1'),ylabel('X2'),zlabel('Y')
    figure(2),clf
    [x11,x33]=meshgrid(x1a:r1/75:x1b,x3a:r3/75:x3b);
    plot3(x1,x3,y,'o')
    stem3(x1,x3,y,'filled')
    if n<250
        text(x1,x3,y+.04*ry,str,'fontsize',12)
    end
    pause(.0001)
    hold on
    y2=nfx(b,x11,m2,x33);
    surf(x11,x33,y2)
    axis([x1a x1b x3a x3b ya yb])
    alpha(.85)
    shading interp
    axis tight
    pause(.0001)
    xlabel('X1'),ylabel('X3'),zlabel('Y')
    figure(3),clf
    plot3(x2,x3,y,'o')
    stem3(x2,x3,y,'filled')
    if n<250
        text(x2,x3,y+.04*ry,str,'fontsize',12)
    end
    pause(.0001)
    hold on
    [x22,x33]=meshgrid(x2a:r2/75:x2b,x3a:r3/75:x3b);
    y3=nfx(b,m1,x22,x33);
    surf(x22,x33,y3)
    axis([x2a x2b x3a x3b ya yb])
    alpha(.85)
    shading interp
    axis tight
    pause(.0001)
    xlabel('X2'),ylabel('X3'),zlabel('Y')
    y1=nfx(b,x1,x2,x3);
    Y=[[1:n]' x1,x2,x3,y,y1];
elseif dim==4
    y1=nfx(b,x1,x2,x3,x4);
    Y=[[1:n]' x1,x2,x3,x4,y,y1];
end
hold off
if dim==1
    y11=nfx(b,x11);y1=nfx(b,x);
    plot(x,y,'ko',x11,y11,'r-','markerfacecolor','k','linewidth',2)
    axis tight
    legend('data','fit','location','best')
elseif dim==2
    y11=nfx(b,x11,x22);
    for i=1:np
        for j=1:nq
            if ~isreal(y11(i,j));
                y11(i,j)=nan;
            end
        end
    end
    y1=nfx(b,x1,x2);
    plot3(x1,x2,y,'o','markerfacecolor','b');
    if n<330
        stem3(x1,x2,y,'filled')
    end
    axis tight
    if n<150
        text(x1,x2,y+.04*ry,str,'fontsize',14)
    end
    hold on
    surf(x11,x22,y11)
    axis tight
    alpha(.85)
    shading interp
    %y1=nfx(b,x1,x2);
    %plot3(x1,x2,y1,'r-','linewidth',3)
    xlabel('X1'),ylabel('X2'),zlabel('Y')
end
db=abs(b)/1e7;
delta=zeros(1,k);
d=y-y1;
for j=1:k
    bl=b;bh=b;
    bl(j)=bl(j)-db(j);bh(j)=bh(j)+db(j);
    if dim==1
        xl=(nfx(bh,x)-nfx(bl,x))/2/db(j);
    elseif dim==2
        xl=(nfx(bh,x1,x2)-nfx(bl,x1,x2))/2/db(j);
    elseif dim==3
        xl=(nfx(bh,x1,x2,x3)-nfx(bl,x1,x2,x3))/2/db(j);
    else
        xl=(nfx(bh,x1,x2,x3,x4)-nfx(bl,x1,x2,x3,x4))/2/db(j);
    end
    xx(:,j)=xl;
    delta(j)=d'*xl;
end
a=xx'*xx;
SEb=sqrt(MSe*diag(inv(a))');
t=bmn./SEb
save out.txt bmn RSS MSe R2 t -ascii