warning off all;
dia=zeros(k);bn=bm;ramda=.013;
for r=1:3000
    b(isnan(b))=bn(isnan(b));
    delta=zeros(1,k);
    if dim==1
        yhat=nfx(b,x);
    elseif dim==2
        yhat=nfx(b,x1,x2);
    elseif dim==3
        yhat=nfx(b,x1,x2,x3);
    else
        yhat=nfx(b,x1,x2,x3,x4);
    end
    d=y-yhat;
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
    if sum(sum(isnan(a)))>=1
        break
    end
    aa=a+ramda*diag(diag(a));
    if det(aa)==0        
        aa=a+.13*diag(diag(a));
        if det(aa)==0
            dia=.8*diag(diag(a));
            aa=a+dia;
            if det(aa)==0
                dia=1.5*diag(diag(a))+.01*rand(k);
                aa=a+dia;
                if det(aa)==0
                    dia=1e3*diag(diag(a))+5*rand(k);
                end
            end
        end
    end
    db=delta/aa;bo=b+.35*db;db=3.5e-3*db;
    if dim==1
        yhat=nfx(bo,x);
    elseif dim==2
        yhat=nfx(bo,x1,x2);
    elseif dim==3
        yhat=nfx(bo,x1,x2,x3);
    else
        yhat=nfx(bo,x1,x2,x3,x4);
    end
    ssr=(y-yhat)'*(y-yhat);    
    if ssr<=qv(v2) & ~isnan(ssr) & isreal(ssr)
        b=bo;bn=bo;
        ramda=max(.1*ramda,eps);
        qv(v2)=ssr;        
        t1=t1+b;s1=s1+b.*b;n1=n1+1;
        if ssr<qm(v1) 
            qm(v1)=ssr;bm=bo;
        end
    elseif isnan(ssr) | ~isreal(ssr)
        b=bn;db=abs(b)/1e7;
        ramda=min(10*ramda,1e8);
        qv(v2)=qm(v1)+c;
    else
        db=.2*db;
        ramda=min(5*ramda,1e3);
        b=b+db;
    end
    v2=v2+1;qv(v2)=qv(v2-1);
    if v2>1000
        for v2=1:10
            qv(v2)=qv(v2+990);
        end
    end
    if v2>12 & (qv(v2-7)-qv(v2))/qv(v2)<=1e-7
        break
    end
end