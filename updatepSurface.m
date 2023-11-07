function [p1,p2,kk]=updatepSurface(p10,p20,s1,s2,tol)

a=0.2;
p1=p10;
p2=p20;

pp1=p1;
pp2=p2;
ss=10;
kk=0;
while ss>tol && kk<50
    kk=kk+1;
    deno=(1+p1.^2+p2.^2).^(5/2);
    fac=s1-s2./deno;
    p1=a*p1+(1-a)*s1*p10./fac;
    
    p2=a*p2+(1-a)*s1*p20./fac;
    
    ss1=max(abs(p1(:)-pp1(:)));
    ss2=max(abs(p2(:)-pp2(:)));
    ss=max(ss1,ss2);
    
    pp1=p1;
    pp2=p2;
end