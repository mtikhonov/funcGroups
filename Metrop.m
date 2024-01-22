function [GpM,Coeff,inE]=Metrop(abd,Y,lin,Nsamples,Nstrains,k) 
% The Metropolis algorithm 
% Input: 
% abd: abundance of each species in each sample.
% Y: function in each sample.
% lin: if lin equals 1, use linear regression, otherwise use quadratic regression. 
% Nsamples: # of samples.
% k: number of groups required for a grouping. (k can be a vector to ask
%    for several groups.)
% 
% Output:
% GpM: the grouping(s)
% Coeff: the coefficients of the regression for each grouping
% inE: in-sample RMSE for each grouping

Nstep=10000;
maxNgroup=min(20,Nstrains-1); % max number of groups
GpM=randGp(Nstrains,1:maxNgroup); % initiate the current list of best groupings with 1,...,maxNgroup groups.
for j=maxNgroup:-1:1
Emin(j,1)=EN(abd,Nsamples,GpM(1,:),Y,lin); % the RMSE of each grouping in the current list
end
recdCef=(nargout==3); % record the coefficients or not
if recdCef
Coeff=cell(length(k),1);
end

for i=1:Nstep   
    candGp=drawEvent(GpM(randsample(maxNgroup,1),:),maxNgroup,2);  % candidate grouping
    [candE,indx,coeff]=EN(abd,Nsamples,candGp,Y,lin);    
    if Emin(indx)>candE %|| rand<exp(-beta*(candE-Emin(indx)))   
        GpM(indx,:)=candGp;                
        Emin(indx)=candE;              
        if recdCef
            idk=indx==k;
            if any(idk)
            Coeff{idk}=coeff;
            end
        end
    end      
end
GpM=GpM(k,:);
if recdCef
    inE=Emin(k);
end
end

function candGp = drawEvent(oriGp,maxNGp,minNGp)  % generate a new candidate grouping.  
    candGp=oriGp;
    NGp=max(oriGp);
    if (NGp<maxNGp && rand>0.5) || NGp<=minNGp % split a group              
        sizeGp=accumarray(candGp',1);
        Gspl=datasample(find(sizeGp>1),1); % group for splitting
        spl=find(candGp==Gspl);      % species in the group   
        Nspl=sizeGp(Gspl);

        g1g2=randperm(Nspl,2);        
        s=rand(1,Nspl)>0.5;
        s(g1g2)=[1 0];   
        candGp(spl(s))=NGp+1;
    else  % combine groups
        rp=randperm(NGp,2);
        candGp(candGp==rp(2))=rp(1);
        candGp(candGp>rp(2))=candGp(candGp>rp(2))-1;
    end 
end

function [E,N,coeff]=EN(abd,Nsamples,Gp,Y,lin) % return the RMSE (E), # of groups (N) and coefficients of regression (coeff) of the grouping Gp. 
N=max(Gp);
cgabd=zeros(Nsamples,N); 
for i = 1:N    
    Gpi=Gp == i;
    cgabd(:,i) = sum(abd(:,Gpi),2);
end
if lin==1
    X=[ones(Nsamples,1),cgabd];
    [coeff,~,r1]=regress(Y,X);
    E=rms(r1);  
else
    D = x2fx(cgabd,'quadratic');   
    [coeff,~,r2]=regress(Y,D);
    E=rms(r2);   
end
end