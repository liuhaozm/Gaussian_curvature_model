clear
noiseLevel=0.01;

f=im2double(imread('peppers.tiff'));
f=imresize(rgb2gray(f),0.5);

[M,N]=size(f);
f0=f;
%%
%%%%%%% Gaussian noise
fn=f0+randn(M,N)*sqrt(noiseLevel);


tol=1e-5;
tau=0.05;
alpha1=0.2;
beta1=0.6;
printerr=1;
[u,tt]=GaussianModel(fn,tau,alpha1,beta1,tol,printerr);

figure
subplot(1,2,1)
imshow(fn)
subplot(1,2,2)
imshow(u)