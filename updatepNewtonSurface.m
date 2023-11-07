function [p1,p2,kk]=updatepNewtonSurface(p10,p20,s1,s2,tol)

dt=1;
p1=p10;
p2=p20;

pp1=p1;
pp2=p2;
ss=10;
kk=0;
while ss>tol && kk<=50
    kk=kk+1;
    deno=(1+p1.^2+p2.^2);
    a11=s1+s2.*(4*p1.^2-p2.^2-1)./deno.^(7/2);
    a12=s2.*5.*p1.*p2./deno.^(7/2);
    a21=a12;
    a22=s1+s2.*(4*p2.^2-p1.^2-1)./deno.^(7/2);
    
    delta=a11.*a22-a12.*a21;
    
    f1=s1.*p1-s2.*p1./deno.^(5/2)-s1*p10;
    f2=s1.*p2-s2.*p2./deno.^(5/2)-s1*p20;
    
    p1=p1-dt*(a22.*f1-a12.*f2)./delta;
    p2=p2-dt*(-a21.*f1+a11.*f2)./delta;
    
    
    ss1=max(abs(p1(:)-pp1(:)));
    ss2=max(abs(p2(:)-pp2(:)));
    ss=max(ss1,ss2);
    
    
    pp1=p1;
    pp2=p2;
end