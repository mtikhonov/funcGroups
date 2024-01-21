function GpE=EQObls2(abd,Y,Nstrains)
opts = optimoptions('ga','Display','none','FunctionTolerance',1e-9,'MaxGenerations',10000,'MaxStallGenerations',500,'PopulationSize',100);
Gp=ga(@(x) obj(x,abd,Y),Nstrains,[],[],[],[],zeros(Nstrains,1),ones(Nstrains,1),[],1:Nstrains,opts);
GpE=Gp+1;
end

function e2 = obj(x,abd,Y)
if ~nnz(x) || all(x)
    A=x2fx(sum(abd,2));
else
    A=x2fx(abd*([x;~x]'));
end
    z=A\Y;
    e=A*z-Y;
    e2=e'*e;
end