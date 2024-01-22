function GpR = randGp(Nstrains,k)
for ii=length(k):-1:1
ball=randperm(Nstrains);
barrier=[0,sort(randperm(Nstrains-1,k(ii)-1))];
gp=k(ii)*ones(1,Nstrains);
for i=1:k(ii)-1
    idx=(barrier(i)+1):barrier(i+1);
    gp(ball(idx))=i;
end
GpR(ii,:)=gp;
end
end