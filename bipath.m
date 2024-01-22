function [abd,Y,expectFine,expectCoarse]=bipath(sizeGp,C,Rr,cSamp,Nsamples,noise)
Ngp=3; idfun=4;
w=[0.4;0.5;0.5;0.4];
[in,Rct2n,Rct2h,sig,m,expectFine,expectCoarse]=layPara(Rr,Ngp,C,sizeGp,w);
Tintrv=[0 100];

b=3; d=1;
K=[4;zeros(3,1)];%3*ones(Rr,1)];
h0=ones(4,1);
for i=Nsamples:-1:1   
n0=zeros(C,1);
n0(randperm(C,cSamp))=1; %exp(randn(10,1));
y0=[n0;h0];
[~,y]=ode45(@(t,y) CrFd(t,y,K,in,Rct2n,Rct2h,w,sig,m,C,b,d),Tintrv,y0);
abd(i,:)=y(end,1:C);
Y(i,1)=y(end,C+idfun);
end
abd=max(abd,0).*max(1+noise*randn(Nsamples,C),0);
Y=max(Y,0).*max(1+noise*randn(Nsamples,1),0);
% eliminate rare species
%sur=max(abd)>0.001*max(abd,[],"all");
%abd=abd(:,sur);
%expectFine=expectFine(sur);
%expectCoarse=expectCoarse(sur);
end

function dy=CrFd(t,y,K,in,Rct2n,Rct2h,w,sig,m,C,b,d)
y=max(y,0);
n=y(1:C);
h=y(C+1:end);
T=sig'*n;
inch=b./(1+T);
dn=n.*(Rct2n*(h(in).*(1-w))+sig*inch-m);
dh=K+Rct2h*((Rct2n'*n).*h(in))-d.*h;
dy=[dn;dh];
end

function [in,Rct2n,Rct2h,sig,m,expectFine,expectCoarse]=layPara(Rr,Ngp,C,sizeGp,w)
expectFine=[repelem(1:4,sizeGp/2),5*ones(1,sizeGp)];
expectCoarse=repelem(1:Ngp,sizeGp);
nfun=2*sizeGp;
in=[1 1 2 3]';%repelem([1:3 0],[10 5 5 80]);
Rct2n=full(sparse(1:nfun,expectFine(1:nfun),1,C,4));
Rct2h=[-1,-1,0,0;w(1),0,-1,0;0,w(2),0,-1;0,0,w(3),w(4)];
sig=rand(C,Rr)<0.3;
for i=find(~sum(sig))
sig(randi(C),i)=1;
end
for i=find(~sum(sig,2))
    sig(i,randi(Rr))=1;
end
m=max(0.01,sum([[repelem(w,sizeGp/2);zeros(sizeGp,1)],sig],2)+0.01*randn);
end