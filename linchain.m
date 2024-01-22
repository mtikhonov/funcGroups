function [abd,Y,expect]=linchain(sizeGp,C,Rr,cSamp,Nstp,idfun,Nsamples,noise)
Ngp=Nstp+1; % # of groups.
[G,sig,D,m,expect]=linPara(Rr,Ngp,sizeGp,C);
Tintrv=[0 100];
b=3; d=1;
K1=2^Nstp;
K=[K1;zeros(Nstp,1)];%3*ones(Rr,1)];
h0=ones(Ngp,1);
for i=Nsamples:-1:1   
n0=zeros(C,1);
n0(randperm(C,cSamp))=1; %exp(randn(10,1));
y0=[n0;h0];

[~,y]=ode45(@(t,y) CrFd(t,y,K,G,sig,D,m,C,b,d),Tintrv,y0);
abd(i,:)=y(end,1:C);
Y(i,1)=y(end,C+idfun);
end
abd=max(abd,0).*max(1+noise*randn(Nsamples,C),0);
Y=max(Y,0).*max(1+noise*randn(Nsamples,1),0);
% eliminate rare species
%sur=max(abd,[],1)>0.001*max(abd,[],"all");
%abd=abd(:,sur);
%expect=expect(sur);
end

function dy=CrFd(t,y,K,G,sig,D,m,C,b,d)
y=max(y,0);
n=y(1:C);
h=y(C+1:end);
T=sig'*n;
inch=b./(1+T);
dn=n.*(((1-sum(D)).*G)*h+sig*inch-m);
dh=K-(G'*n).*h+D*((G'*n).*h)-d.*h;
dy=[dn;dh];
end

function [G,sig,D,m,expect]=linPara(Rr,Ngp,sizeGp,C)
expect=repelem(1:Ngp,sizeGp);
nfun=C-sizeGp(end);
G=full(sparse(1:nfun,expect(1:nfun),1,C,Ngp));
sig=rand(C,Rr)<0.3;
for i=find(~sum(sig))
sig(randi(C),i)=1;
end
for i=find(~sum(sig,2))
    sig(i,randi(Rr))=1;
end
w=0.5;
D=diag(w*ones(1,Ngp-1),-1);
m=max(0.01,sum([w*G,sig],2)+0.01*randn);
end