% MATLAB code for "An Operator-Splitting Method for the Gaussian Curvature 
% Regularization Model with Applications to Surface Smoothing and Imaging" 
% by Liu, Hao and Tai, Xue-Cheng and Glowinski, Roland.
% https://arxiv.org/abs/2108.01914
%
% If this code is useful, pleae cite our paper.
%

% Copyright (c) 2022 Hao Liu (haoliu AT hkbu.edu.hk)
% Department of Mathematics,
% Hong Kong Baptist University
% Kowloon, Hong Kong
% https://www.math.hkbu.edu.hk/~haoliu/
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.


function [u,tt]=GaussianModel(fn,tau,alpha1,beta1,tol,printerr)
% Arguments:
% fn: noisy image
% tau: time step
% alpha1: weight of the TV term
% beta1: (1/2beta1) is the weight of the fidelity term
% tol: stopping critetion
% printerr: whether print error and energy history.
%%
[M,N]=size(fn);
h=1;

xx=(1:M)-1;
yy=(1:N)-1;
[y,x]=meshgrid(yy,xx);
zx=2*pi*x/M;
zy=2*pi*y/N;

gamma=1;


w0=((1-exp(-sqrt(-1)*zx)).*(exp(sqrt(-1)*zx)-1)+...
    (1-exp(-sqrt(-1)*zy)).*(exp(sqrt(-1)*zy)-1));


u=fn;


[p1,p2]=gradu(u,h,1);


[H11,H12]=gradu(p1,h,2);
[H21,H22]=gradu(p2,h,2);


detH=H11.*H22-H12.*H21;


kk=0;
ss=10;

tic
while ss>tol
    kk=kk+1;
    
    s1=gamma;
    s2=3*tau*(abs(H11.*H22-H12.*H21));
    
    [p1_1,p2_1,countp]=updatepSurface(p1,p2,s1,s2,1e-5);
%     [p1_1,p2_1,countp]=updatepNewtonSurface(p1,p2,s1,s2,1e-5);

    s1=tau./(1+p1_1.^2+p2_1.^2).^(3/2);
    
    [H11_1,H22_1,H12_1,H21_1,countH]=updateHAlter(H11,H22,H12,H21,s1,1e-5);
    
    p1_1=max(0,1-tau*alpha1./(sqrt(p1_1.^2+p2_1.^2+1e-10))).*p1_1;
    p2_1=max(0,1-tau*alpha1./(sqrt(p1_1.^2+p2_1.^2+1e-10))).*p2_1;
    
    %%% step 2
    
    [H11x_1,H11y_1]=gradu(H11_1,h,1);
    [H22x_1,H22y_1]=gradu(H22_1,h,1);
    [H12x_1,H12y_1]=gradu(H12_1,h,1);
    [H21x_1,H21y_1]=gradu(H21_1,h,1);
    
    
    divH1=H11x_1+H12y_1;
    divH2=H21x_1+H22y_1;
    
    rhs1=(gamma*p1_1-divH1)*h^2;
    rhs2=(gamma*p2_1-divH2)*h^2;
    rhs1F=fft2(rhs1);
    rhs2F=fft2(rhs2);
    
    p1_2=real(ifft2((rhs1F)./(gamma*h^2-w0)));
    p2_2=real(ifft2((rhs2F)./(gamma*h^2-w0)));
    
    
    [H11_2,H12_2]=gradu(p1_2,h,2);
    [H21_2,H22_2]=gradu(p2_2,h,2);
    
    %%%% step 3
    [p1x_2,p1y_2]=gradu(p1_2,h,2);
    [p2x_2,p2y_2]=gradu(p2_2,h,2);
    
    divp=p1x_2+p2y_2;
    rhs=h^2*divp*gamma-tau*h^2*fn/beta1;
    rhsF=fft2(rhs);
    
    w=w0*h^2*(gamma)-tau*h^2/beta1;
    u_3=real(ifft2(rhsF./w));
    [p1_3,p2_3]=gradu(u_3,h,1);
    
    H11_3=H11_2;
    H12_3=H12_2;
    H21_3=H21_2;
    H22_3=H22_2;
    
    %%%%%%% compute energy
    uu=u_3;
    
    [pp1,pp2]=gradu(uu,h,1);
    
    [HH11,HH12]=gradu(pp1,h,2);
    [HH21,HH22]=gradu(pp2,h,2);
    
    HH12ex=expandf(HH12);
    HH12=(HH12ex(2:end-1,2:end-1)+HH12ex(1:end-2,3:end) )/2;
    HH21ex=expandf(HH21);
    HH21=(HH21ex(2:end-1,2:end-1)+HH21ex(3:end,1:end-2))/2;
    
    detHH=HH11.*HH22-HH12.*HH21;
    
    ener1=sum(abs(detHH)./((1+pp1.^2+pp2.^2).^2),'all');
    ener2=1/2*sum((fn-uu).^2,'all')/beta1;
    ener3=alpha1*(sum(sqrt(pp1.^2+pp2.^2),'all'));
    ener=ener1+ener2+ener3;
    
    ss=sqrt(sum((u(:)-u_3(:)).^2))/sqrt(sum((u_3(:)).^2));
    if printerr
        sprintf('%i %2i %2i %e %e',kk,countp, countH,ener,ss)
    end
    
    %%%% update
    u=u_3;
    p1=p1_3;
    p2=p2_3;
    H11=H11_3;
    H12=H12_3;
    H21=H21_3;
    H22=H22_3;
    
end

tt=toc;


