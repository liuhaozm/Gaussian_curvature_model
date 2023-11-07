function pexpand=expandf(p)

[M,N]=size(p);
pexpand=zeros(M+2,N+2);
pexpand(2:end-1,2:end-1)=p;
pexpand(1,:)=pexpand(end-1,:);
pexpand(end,:)=pexpand(2,:);
pexpand(:,1)=pexpand(:,end-1);
pexpand(:,end)=pexpand(:,2);
