function [abd,Y,expect]=branch(sizeGp,C,Rr,cSamp,Nsamples,noise)
Ngp=4; idfun=3; w=0.5;
[G,sig,D,m,expect]=brchPara(Rr,Ngp,sizeGp,C,w);
Tintrv=[0 100];

b=3; d=1;
K=[4;zeros(3,1)];%3*ones(Rr,1)];
h0=ones(Ngp,1);
for i=Nsamples:-1:1   
n0=zeros(C,1);
n0(randperm(C,cSamp))=1; %exp(randn(10,1));
y0=[n0;h0];
[~,y]=ode45(@(t,y) CrFd(t,y,K,G,sig,D,m,C,b,w,d),Tintrv,y0);
abd(i,:)=y(end,1:C);
Y(i,1)=y(end,C+idfun);
end
abd=max(abd,0).*max(1+noise*randn(Nsamples,C),0);
Y=max(Y,0).*max(1+noise*randn(Nsamples,1),0);
% eliminate rare species
%sur=max(abd)>0.001*max(abd,[],"all");
%abd=abd(:,sur);
%expect=expect(sur);
end

function dy=CrFd(t,y,K,G,sig,D,m,C,b,w,d)
y=max(y,0);
n=y(1:C);
h=y(C+1:end);
T=sig'*n;
inch=b./(1+T); in=[1;2;2];
dn=n.*(((1-w).*G)*h(in)+sig*inch-m);
dh=K+D*((G'*n).*h(in))-d.*h;
dy=[dn;dh];
end

function [G,sig,D,m,expect]=brchPara(Rr,Ngp,sizeGp,C,w)
expect=repelem(1:Ngp,sizeGp);
nfun=C-sizeGp(end);
G=full(sparse(1:nfun,expect(1:nfun),1,C,Ngp-1));
sig=rand(C,Rr)<0.3;
for i=find(~sum(sig))
sig(randi(C),i)=1;
end
for i=find(~sum(sig,2))
    sig(i,randi(Rr))=1;
end
D=[-1,0,0;w,-1,-1;0,w,0;0,0,w];
m=max(0.01,sum([(1-w)*G,sig],2)+0.01*randn);
end