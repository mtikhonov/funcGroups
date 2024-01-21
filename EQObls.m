function GpE=EQObls(abd,Y,Nstrains,Nsamples)
opts = optimoptions('ga','Display','none','MaxGenerations',10000, ...
    'FunctionTolerance',1e-9,'MaxStallGenerations',500, ...
    'PopulationSize',100);
Gp=ga(@(x) obj(x,abd,Y,Nsamples),Nstrains,[],[],[],[],zeros(Nstrains,1),ones(Nstrains,1),[],1:Nstrains,opts);
GpE=Gp+1;
end

function aic = obj(x,abd,Y,Nsamples)
if any(x)
    A=[ones(Nsamples,1),abd*x'];
    [~,~,e]=regress(Y,A);
    e2=e'*e;
    aic=Nsamples*log(e2/Nsamples)+2*nnz(x);
else
    aic=Nsamples*log(var(Y,1));
end
end

