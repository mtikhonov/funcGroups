function score=JacSim(expect,Gp,k)
% Calculating the recovery quality of true groups indicated by "k".
% If "k" is not a input, calculate for all true groups.
if nargin==2
k=unique(expect);
end
x=numel(k);
for ii=size(Gp,1):-1:1
    gp=Gp(ii,:);
    y=max(gp);
    sc=zeros(y,x);
    for i=1:x
        for j=1:y
            a=(expect==k(i));
            b=(gp==j);
            sc(j,i)=jaccard(a,b);
        end
    end
    score(ii,:)=max(sc);
end
end