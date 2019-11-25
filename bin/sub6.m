b(isnan(b))=bmn(isnan(b));
if ~isreal(b)==1,b=bmn;end
if dim==1
    yhat=nfx(b,x);
elseif dim==2
    yhat=nfx(b,x1,x2);
elseif dim==3
    yhat=nfx(b,x1,x2,x3);
else
    yhat=nfx(b,x1,x2,x3,x4);
end
yhat(isnan(yhat))=y(isnan(yhat))*1000;
q=(y-yhat)'*(y-yhat);
if q-(.9+.18*v)*c<=q1
    t1=t1+b;s1=s1+b.*b;n1=n1+1;nc=nc+1;q1=q;bo=b;
    if q<=qm(v1)
        bm=b;qm(v1)=q;b0=b;
    end
end