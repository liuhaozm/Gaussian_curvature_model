function [H11,H22,H12,H21,count]=updateHAlter(H110,H220,H120,H210,s1,tol)

fac=0.8;
eps2=1e-14;
H11=H110;
H12=H120;
H22=H220;
H21=H210;

HH11=H11;
HH12=H12;
HH21=H21;
HH22=H22;

sss=10;
count=0;
while sss>tol && count<50
    count=count+1;
a=H22;
b=H21;
alpha1=H110;
beta1=H120;
cri1=(a.*alpha1-b.*beta1)-(a.^2.*s1+b.^2.*s1);
cri2=(a.*alpha1-b.*beta1)+(a.^2.*s1+b.^2.*s1);
%%
HH11=(a==0).*(alpha1)+(a~=0).*(b==0).*(...
    max(0,1-s1.*abs(a)./(abs(alpha1)+1e-8)).*alpha1...
    )+(a~=0).*(b~=0).*(...
    (alpha1-s1.*a).*(cri1>0)+(alpha1+s1.*a).*(cri2<0)+...
    (alpha1.*b.^2+a.*beta1.*b)./(b.^2+(a).^2+eps2).*(cri1<=0).*(cri2>=0)...
    );

HH12=(b==0).*(beta1)+(b~=0).*(a==0).*(...
    max(0,1-s1.*abs(b)./(abs(beta1)+1e-8)).*beta1...
    )+(b~=0).*(a~=0).*(...
    (beta1+s1.*b).*(cri1>0)+(beta1-s1.*b).*(cri2<0)+...
    (beta1.*a.^2+b.*alpha1.*a)./(a.^2+(b).^2+eps2).*(cri1<=0).*(cri2>=0)...
    );

a=fac*HH11+(1-fac)*H11;
b=fac*HH12+(1-fac)*H12;
alpha1=H220;
beta1=H210;
cri1=(a.*alpha1-b.*beta1)-(a.^2.*s1+b.^2.*s1);
cri2=(a.*alpha1-b.*beta1)+(a.^2.*s1+b.^2.*s1);
HH22=(a==0).*(alpha1)+(a~=0).*(b==0).*(...
    max(0,1-s1.*abs(a)./(abs(alpha1)+1e-8)).*alpha1...
    )+(a~=0).*(b~=0).*(...
    (alpha1-s1.*a).*(cri1>0)+(alpha1+s1.*a).*(cri2<0)+...
    (alpha1.*b.^2+a.*beta1.*b)./(b.^2+(a).^2+eps2).*(cri1<=0).*(cri2>=0)...
    );

HH21=(b==0).*(beta1)+(b~=0).*(a==0).*(...
    max(0,1-s1.*abs(b)./(abs(beta1)+1e-8)).*beta1...
    )+(b~=0).*(a~=0).*(...
    (beta1+s1.*b).*(cri1>0)+(beta1-s1.*b).*(cri2<0)+...
    (beta1.*a.^2+b.*alpha1.*a)./(a.^2+(b).^2+eps2).*(cri1<=0).*(cri2>=0)...
    );


ss1=max(abs(H11(:)-HH11(:)));
ss2=max(abs(H22(:)-HH22(:)));
ss3=max(abs(H12(:)-HH12(:)));
ss4=max(abs(H21(:)-HH21(:)));
sss=max([ss1,ss2,ss3,ss4]);

H11=fac*HH11+(1-fac)*H11;
H22=fac*HH22+(1-fac)*H22;
H12=fac*HH12+(1-fac)*H12;
H21=fac*HH21+(1-fac)*H21;
end