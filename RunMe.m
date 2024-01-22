%% 
% Codes for generating Fig. 2-4 and Fig. S2-5 of
% 
% Y. Zhao, O. X. Cordero and M. Tikhonov, Linear-regression-based algorithms 
% can succeed at identifying microbial functional groups despite the nonlinearity 
% of ecological function.

clear
Nsamples=900; noise=0.1;
Nsamples_array=100:100:900;
noise_array=0.05:0.05:0.3;
Ntest=50; % Number of datasets generated for testing. One can change it to a smaller value for faster running.
MakeFig2_S4_S5(Ntest,Nsamples,noise,Nsamples_array,noise_array) % validation of 3 algorithms for N=3
MakeFig3(Ntest,Nsamples,noise) % longer chain
MakeFig4(2*Ntest,Nsamples_array,noise_array) % linear vs. quadratic 
MakeFigS2(Ntest,Nsamples,noise) % EQO-2g
MakeFigS3(Ntest,Nsamples,noise) % other functions
%%
function MakeFig2_S4_S5(Ntest,Nsamples,noise,Nsamples_array,noise_array)
Nstp=2; idfun=3;
sizeGp=16; C=48;
cSamp=15; Rr=15;
k=2:4;
if exist("fig2data.mat",'file')
    load("fig2data.mat")
else
fprintf('Generating Fig 2 data.\n')
% Panel A & B
fprintf('Panels A, B.\n')
parfor ii=1:Ntest   
    [abd,Y,expect]=linchain(sizeGp,C,Rr,cSamp,Nstp,idfun,Nsamples,noise);    
    Nstrains=size(abd,2);
    gpK=KMeans(abd,Y,k); GpK(ii,:,:)=gpK';
    gpM=Metrop(abd,Y,1,Nsamples,Nstrains,k); GpM(ii,:,:)=gpM';              
    gpE=EQObls(abd,Y,Nstrains,Nsamples); GpE(ii,:)=gpE;
    gpR=randGp(Nstrains,k);
    Sc(ii,:,:)=JacSim(expect,[gpE;gpK;gpM;gpR]);
end

% Panel C
fprintf('Panel C.\n')
NS=length(Nsamples_array); NO=length(noise_array);
parfor ii=1:Ntest
    for ns=1:NS  
      for no=1:NO  
        [abd,Y,expect]=linchain(sizeGp,C,Rr,cSamp,Nstp,idfun,Nsamples_array(ns),noise_array(no));    
        Nstrains=size(abd,2); %Expect{ns,no,ii}=expect;
        gpK=KMeans(abd,Y,k); %GpK{ns,no,ii}=gpK;   
        ScoreK(ns,no,:,ii)=mean(JacSim(expect,gpK),2); 
        gpM=Metrop(abd,Y,1,Nsamples_array(ns),Nstrains,k); %GpM{ns,no,ii}=gpM;              
        ScoreM(ns,no,:,ii)=mean(JacSim(expect,gpM),2);
        gpE=EQObls(abd,Y,Nstrains,Nsamples_array(ns)); %GpE{ns,no,ii}=gpE;          
        ScoreE(ns,no,:,ii)=mean(JacSim(expect,gpE),2);   
       end              
    end    
end
fprintf('Saving Fig 2 data.\n')
save fig2data Sc Nsamples_array noise_array GpK ScoreK GpM ScoreM GpE ScoreE
fprintf('Done.\n') 
end
plotFig2(Sc,2/3,ScoreM,Nsamples_array,noise_array)
plotFigS4(GpE,GpK,GpM)
plotFigS5(ScoreE,ScoreK,ScoreM,noise,Nsamples)
end

function MakeFig3(Ntest,Nsamples,noise)
Nstp=3; idfun=4;
sizeGp=12; C=48;
cSamp=15; Rr=15; k=2:5;
if exist("fig3Adata.mat",'file')
    load("fig3Adata.mat")
else
    fprintf('Generating Fig 3A data.\n')
    parfor ii=1:Ntest   
        [abd,Y,expect]=linchain(sizeGp,C,Rr,cSamp,Nstp,idfun,Nsamples,noise);    
        Nstrains=size(abd,2);
        gpM=Metrop(abd,Y,1,Nsamples,Nstrains,k); 
        gpR=randGp(Nstrains,k);
        ScM(:,:,ii)=JacSim(expect,gpM,1:4);
        ScR(:,:,ii)=JacSim(expect,gpR,1);
    end
    fprintf('Saving Fig 3A data.\n')
    save fig3Adata ScM ScR k
    fprintf('Done.\n')
end

sizeGp=[24,16,12,10,8]; C=[48,48,48,50,48];
if exist("fig3Bdata.mat",'file')
    load("fig3Bdata.mat")
else
    fprintf('Generating Fig 3B data.\n')
    parfor ii=1:Ntest
        fprintf('%d/%d\n',ii,Ntest);
        for stp=1:5
            idfun=stp+1;  %sizeGp=C/idfun;
            [abd,Y,expect]=linchain(sizeGp(stp),C(stp),Rr,cSamp,stp,idfun,Nsamples,noise);                  
            Nstrains=size(abd,2);
            % GpK=KMeans(abd,Y,[2 idfun]); farK nearK   
            GpM=Metrop(abd,Y,1,Nsamples,Nstrains,[2 idfun]);          
            GpR=randGp(Nstrains,[2 idfun]);       
            js=JacSim(expect,[GpM;GpR]);  
            % farK(ii,stp)=sc(2,1); nearK(ii,stp)=sc(1,2);
            ScoreMR(stp,:,ii)=mean(js([2,4],:),2);
            farM(ii,stp)=js(2,1); nearM(ii,stp)=js(1,stp);   
            farR(ii,stp)=js(4,1); nearR(ii,stp)=js(3,stp);
        end      
    end
    fprintf('Saving Fig 3B data.\n')
    save fig3Bdata ScoreMR farM nearM farR nearR Ntest
    fprintf('Done.\n')
end
plotFig3(ScM,ScR,k,Ntest,nearR,farR,nearM,farM)
end

function MakeFig4(Ntest,Nsamples_array,noise_array)
Nstp=2; idfun=3;
sizeGp=16; C=48;
cSamp=15; Rr=15;
if exist("fig4data.mat",'file')
    load("fig4data.mat")
else
fprintf('Generating Fig 4 data.\n')
NS=length(Nsamples_array); NO=length(noise_array);
parfor ii=1:Ntest 
    fprintf('%d/%d.\n',ii,Ntest);      
    for ns=1:NS 
      for no=1:NO           
        [abd,Y,expect]=linchain(sizeGp,C,Rr,cSamp,Nstp,idfun,2*Nsamples_array(ns),noise_array(no));    
        Nstrains=size(abd,2);
        abd1=abd(1:Nsamples_array(ns),:); Y1=Y(1:Nsamples_array(ns),:);
        abd2=abd(Nsamples_array(ns)+1:end,:); Y2=Y(Nsamples_array(ns)+1:end,:); mY2=mean(Y2);      
        
        [gpl,CoeffL,~]=Metrop(abd1,Y1,1,Nsamples_array(ns),Nstrains,3);    
        lin(ns,no,ii)=mean(JacSim(expect,gpl)); GpMl{ns,no,ii}=gpl;
        cgabdl=(groupsummary(abd2',gpl(end,:)','sum'))'; Yprl=x2fx(cgabdl)*CoeffL{end};
        R2lin(ns,no,ii)=1-sum((Yprl-Y2).^2)/sum((Y2-mY2).^2);
        
        warning('off','stats:regress:RankDefDesignMat')
        [gpq,CoeffQ,~]=Metrop(abd1,Y1,2,Nsamples_array(ns),Nstrains,3);        
        qua(ns,no,ii)=mean(JacSim(expect,gpq)); GpMq{ns,no,ii}=gpq; 
        cgabdq=(groupsummary(abd2',gpq(end,:)','sum'))'; Yprq=x2fx(cgabdq,"quadratic")*CoeffQ{end};
        R2qua(ns,no,ii)=1-sum((Yprq-Y2).^2)/sum((Y2-mY2).^2);       
        warning('on','stats:regress:RankDefDesignMat')
       end
    end 
end
fprintf('Saving Fig 4 data.\n')
save fig4data lin qua R2qua R2lin GpMl GpMq Nsamples_array noise_array   
fprintf('Done.\n')
end
plotFig4(lin,qua,R2lin,R2qua,Nsamples_array,noise_array)
end

function MakeFigS2(Ntest,Nsamples,noise)
Nstp=1; idfun=2;
sizeGp=24; C=48;
cSamp=15; Rr=15; k=2;
if exist("figS2data.mat",'file')
    load("figS2data.mat")
else
    fprintf('Generating Fig S2 data.\n')
    [abd1,Y1,expect1]=linchain(sizeGp,C,Rr,cSamp,Nstp,idfun,Nsamples,noise);
    parfor ii=1:Ntest   
        [abd,Y,expect]=linchain(sizeGp,C,Rr,cSamp,Nstp,idfun,Nsamples,noise); 
        Nstrains=size(abd,2);
        gpK=KMeans(abd,Y,k);
        gpM=Metrop(abd,Y,1,Nsamples,Nstrains,k);                
        gpE2=EQObls2(abd,Y,Nstrains);
        gpE=EQObls(abd,Y,Nstrains,Nsamples);
        gpR=randGp(Nstrains,k);
        Sc1(ii,:,:)=JacSim(expect,[gpE;gpE2;gpK;gpM;gpR]);
    end
    parfor ii=1:Ntest   
        [abd,Y,expect]=linchain(sizeGp,C,2*Rr,cSamp,Nstp,idfun,Nsamples,noise); 
        Nstrains=size(abd,2);
        gpK=KMeans(abd,Y,k);
        gpM=Metrop(abd,Y,1,Nsamples,Nstrains,k);                
        gpE2=EQObls2(abd,Y,Nstrains);
        gpE=EQObls(abd,Y,Nstrains,Nsamples);
        gpR=randGp(Nstrains,k);
        Sc2(ii,:,:)=JacSim(expect,[gpE;gpE2;gpK;gpM;gpR]);
    end
    fprintf('Saving Fig S2 data.\n')
    save figS2data abd1 Y1 expect1 Sc1 Sc2
    fprintf('Done.\n')
end
plotFigS2(abd1,Y1,expect1,Sc1,Sc2)
end

function MakeFigS3(Ntest,Nsamples,noise)
sizeGp=16; C=48;
cSamp=15; Rr=15; k=2:5;
if exist("figS3m1data.mat",'file')
    load("figS3m1data.mat")
else
    fprintf('Generating Fig S3m1 data.\n')
    parfor ii=1:Ntest      
        [abd,Y,expect]=linchain(sizeGp,C,Rr,cSamp,2,2,Nsamples,noise);    
        Nstrains=size(abd,2);
        gpM=Metrop(abd,Y,1,Nsamples,Nstrains,k);                
        gpR=randGp(Nstrains,k);
        ScInM(:,:,ii)=JacSim(expect,gpM,1:2);
        ScInR(:,:,ii)=JacSim(expect,gpR,1);
    end
    save figS3m1data ScInM ScInR k
    fprintf('Done.\n')
end
plotFig5(cat(2,ScInM,ScInR),k,{'b^-','ro-','k*-'},Ntest,"figS3A.svg")

sizeGp=12; C=48;
cSamp=15; Rr=15; k=2:5;
if exist("figS3m2data.mat",'file')
    load("figS3m2data.mat")
else
    fprintf('Generating Fig S3m2 data.\n')
    parfor ii=1:Ntest           
        [abd,Y,expect]=branch(sizeGp,C,Rr,cSamp,Nsamples,noise);    
        Nstrains=size(abd,2);
        gpM=Metrop(abd,Y,1,Nsamples,Nstrains,k);                
        gpR=randGp(Nstrains,k);
        ScBrM(:,:,ii)=JacSim(expect,gpM,1:3);
        ScBrR(:,:,ii)=JacSim(expect,gpR,1);       
    end
    save figS3m2data ScBrM ScBrR k
    fprintf('Done.\n')
end
plotFig5(cat(2,ScBrM,ScBrR),k,{'b^-','ro-','gx-','k*-'},Ntest,"figS3B.svg")
 
sizeGp=16; C=48;
cSamp=15; Rr=15; k=2:5;
if exist("figS3m3data.mat",'file')
    load("figS3m3data.mat")
else
    fprintf('Generating Fig S3m3 data.\n')
    parfor ii=1:Ntest       
        [abd,Y,expectFine,expectCoarse]=bipath(sizeGp,C,Rr,cSamp,Nsamples,noise);    
        Nstrains=size(abd,2);    
        gpM=Metrop(abd,Y,1,Nsamples,Nstrains,k);                
        gpR=randGp(Nstrains,k);
        ScBpMc(:,:,ii)=JacSim(expectCoarse,gpM,1:2);
        ScBpMf(:,:,ii)=JacSim(expectFine,gpM,1:4);
        ScBpR(:,:,ii)=JacSim(expectFine,gpR,[1 5]);      
    end
    save figS3m3data ScBpMc ScBpMf ScBpR k
    fprintf('Done.\n')
end
plotFig5(cat(2,ScBpMc,ScBpR(:,2,:)),k,{'b^-','ro-','k*-'},Ntest,"figS3C1.svg")
plotFig5(cat(2,ScBpMf,ScBpR(:,1,:)),k,{'b^--','b^-.','ro--','ro-.','k*-'},Ntest,"figS3C2.svg")
end

function plotFig2(Sc,maxSc2,ScoreM,Nsamples,noise)
W = 50; H = 20;
Sc=mean(Sc,3);
clf
f2=tiledlayout(1,3,"Padding","compact");
set(gcf,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperSize',[W H],'PaperPosition',[0 0 W H], ...
    'Units','Centimeters','Position',[2 2 W H]);
nexttile
hold on
boxchart(Sc(:,[1 2 5 8]),'MarkerSize',10,'MarkerColor',[0.3,0.3,0.3],'MarkerStyle','.','JitterOutliers','on','BoxEdgeColor',"k","BoxFaceColor",[0.5 0.5 0.5])
plot(1:4,mean(Sc(:,[1 2 5 8])),'pentagramk','MarkerSize',12)
hold off
yline(maxSc2,'k--',"LineWidth",1.5)
xticklabels(["EQO","K-Means","Metropolis","Random"])
ylabel("Grouping score")
title("2-group")
ylim([0.2 1])
set(gca,"FontSize",20)
box on
axis square

nexttile
hold on
boxchart(Sc(:,[3 6 9]),'MarkerSize',10,'MarkerColor',[0.3,0.3,0.3],'MarkerStyle','.','JitterOutliers','on','BoxEdgeColor',"k","BoxFaceColor",[0.5 0.5 0.5])
plot(1:3,mean(Sc(:,[3 6 9])),'pentagramk','MarkerSize',12)
yline(maxSc2,'k--',"LineWidth",1.5)
xticklabels(["K-Means","Metropolis","Random"])
xtickangle(30)
ylabel("Grouping score")
title("3-group")
box on
axis square
set(gca,"FontSize",20)

nexttile
me=mean(ScoreM,4);
me3=me(:,:,2);
imagesc(me3)
text(2,size(me3,1),'*','FontSize',35,'HorizontalAlignment','center','VerticalAlignment','middle')
yticks(1:length(Nsamples))
yticklabels(Nsamples)
ylabel('Number of samples')
xticks(1:length(noise))
xticklabels(noise)
xlabel('Relative noise')
title('Score of 3-group Metropolis')
colorbar
axis xy
axis square
set(gca,"FontSize",20)

annotation("textbox",[0.01 0.93 0.01 0.01],'String',"A","LineStyle",'none','FontSize',30)
annotation("textbox",[0.33 0.93 0.01 0.01],'String',"B","LineStyle",'none','FontSize',30)
annotation("textbox",[0.65 0.93 0.01 0.01],'String',"C","LineStyle",'none','FontSize',30)

exportgraphics(f2,'fig2.pdf',"ContentType","vector")
end

function plotFig3(ScM,ScR,k,Ntest,nearR,farR,nearM,farM)
W = 33;
H = 20;
clf
tiledlayout(1,2,"Padding","compact")
set(gcf,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperSize',[W H],'PaperPosition',[0 0 W H], ...
    'Units','Centimeters','Position',[2 2 W H]);
nexttile
plotFig5(cat(2,ScM(:,1:3,:),ScR),k,{'bo--','go-.','ro-','k^-'},Ntest)
ylabel("Recovery quality")
ylim([0.1 1])
text(2.1,0.92,"Direct producers (#3)",'Color','r',"FontSize",18)
text(3.05,0.73,"More upstream (#2)",'Color','g',"FontSize",18)
text(3.5,0.47,"Most upstream (#1)",'Color','b',"FontSize",18)
box on
axis square
set(gca,"FontSize",20)

nexttile
x=2:size(nearM,2)+1;
xconf=[x,x(end:-1:1)];   
y2=mean(nearM);
sem2=std(nearM)/sqrt(Ntest);
yconf2=[y2+sem2,y2(end:-1:1)-sem2(end:-1:1)];
y4=mean(farM);
sem4=std(farM)/sqrt(Ntest);
yconf4=[y4+sem4,y4(end:-1:1)-sem4(end:-1:1)];
y5=mean(farR);
sem5=std(farR)/sqrt(Ntest);
yconf5=[y5+sem5,y5(end:-1:1)-sem5(end:-1:1)];
y6=mean(nearR);
sem6=std(nearR)/sqrt(Ntest);
yconf6=[y6+sem6,y6(end:-1:1)-sem6(end:-1:1)];
hold on
plot(x,y2,'ro-','LineWidth',1.2);
plot(x,y4,'bo--','LineWidth',1.2);
plot(x,y6,'k^-','LineWidth',1.2);
plot(x,y5,'k^--','LineWidth',1.2);
fill(xconf,yconf2,'r',"EdgeColor","none","FaceAlpha",0.2)
fill(xconf,yconf4,'b',"EdgeColor","none","FaceAlpha",0.2)
fill(xconf,yconf5,'k',"EdgeColor","none","FaceAlpha",0.2)
fill(xconf,yconf6,'k',"EdgeColor","none","FaceAlpha",0.2)
hold off
xlabel('Length of chain (N)')
ylabel("Recovery quality")
ylim([0.1 1])
text(4,0.92,"Direct producers (#N-1)",'Color','r','Rotation',-42,"FontSize",18)
text(2.8,0.77,"Most upstream (#1)",'Color','b','Rotation',-42,"FontSize",18)
axis square
yticks(0.2:0.2:1)
xticks(x)
box on
set(gca,"FontSize",20)
saveas(gcf,"fig3.svg")
end

function plotFig4(lin,qua,R2lin,R2qua,Nsamples,noise)
b=[(0:0.01:0.99)',(0:0.01:0.99)',ones(100,1)];
r=[ones(100,1),(0.99:-0.01:0)',(0.99:-0.01:0)'];
bwr=[b;ones(1,3);r];

difS=mean(lin,3)-mean(qua,3);
mxS=min(0.3,max(abs(difS),[],'all'));

difE=mean(R2lin,3)-mean(R2qua,3);
mxE=min(0.3,max(abs(difE),[],'all'));

W = 33;
H = 20;
clf
tiledlayout(1,2,"Padding","compact")
set(gcf,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperSize',[W H],'PaperPosition',[0 0 W H], ...
    'Units','Centimeters','Position',[2 2 W H]);
nexttile
imagesc(difE)
clim([-mxE mxE])
colormap(bwr)
yticks(1:length(Nsamples))
yticklabels(Nsamples)
ylabel('Number of samples')
xticks(1:length(noise))
xticklabels(noise)
xlabel('Relative noise')
if mxE==0.3
colorbar("Ticks",(-0.3):0.1:0.3,'TickLabels',["\leq -0.3",-0.2:0.1:0.2,"\geq 0.3"])
end
% Ploting the phase boundaries. The exact coordinates should be adjust each time manually. 
% hold on
% plot([0.5,1.5,1.5,3.5,3.5],[3.5,3.5,8.5,8.5,9.5],'k--')
% plot([0.5 1.5 1.5 2.5 2.5 3.5 3.5 4.5 4.5 5.5 5.5],[2.5 2.5 3.5 3.5 4.5 4.5 6.5 6.5 8.5 8.5 9.5],'k--')
% hold off
title(["      Ability to predict function","(\Delta out-of-sample R^2)"],"HorizontalAlignment","center")
axis xy
axis square
set(gca,"FontSize",20)

nexttile
imagesc(difS)
clim([-mxS mxS])
colormap(gca,bwr)
if mxS==0.3
colorbar("Ticks",(-0.3):0.1:0.3,'TickLabels',["\leq -0.3",-0.2:0.1:0.2,"\geq 0.3"])
end
yticks(1:length(Nsamples))
yticklabels(Nsamples)
ylabel('Number of samples')
xticks(1:length(noise))
xticklabels(noise)
xlabel('Relative noise')
colorbar
title(["    Ability to detect groups","    (\Delta grouping score)"],"HorizontalAlignment","center")
axis xy
axis square
set(gca,"FontSize",20);
hold on
scatter([],[],20,"red","filled","s")
scatter([],[],10,"blue","filled","s")
% Same as above.
% plot([0.5,1.5,1.5,3.5,3.5],[3.5,3.5,8.5,8.5,9.5],'k--')
% plot([0.5 1.5 1.5 2.5 2.5 3.5 3.5 4.5 4.5 5.5 5.5],[2.5 2.5 3.5 3.5 4.5 4.5 6.5 6.5 8.5 8.5 9.5],'k--')
hold off
[lg,ob]=legend(["Linear is better","Quadratic is better"],"Orientation","horizontal");
lg.Layout.Tile='north';
obj=findobj(ob,'type','patch');
set(obj,"markersize",20)

annotation("textbox",[0.01 0.83 0.01 0.01],'String',"A","LineStyle",'none','FontSize',30)
annotation("textbox",[0.52 0.83 0.01 0.01],'String',"B","LineStyle",'none','FontSize',30)
exportgraphics(gcf,'fig4.pdf',"ContentType","vector")
end

function plotFig5(Sc,x,mark,Ntest,name) % Function for ploting Fig. 3A and all panels of Fig. S3.
xconf=[x,x(end:-1:1)]; 
if nargin==5
figure
end
hold on
for i=1:length(mark)
    data=Sc(:,i,:);
    y=(mean(data,3))';
    sem=(std(data,0,3)/sqrt(Ntest))';
    yconf=[y+sem,y(end:-1:1)-sem(end:-1:1)];
    plot(x,y,mark{i},'LineWidth',1.2);
    fill(xconf,yconf,mark{i}(1),"EdgeColor","none","FaceAlpha",0.2)
end
hold off
xlabel("# of groups in output (k)")
ylabel('Recovery quality')
xticks(x)
set(gca,"FontSize",20)
axis square
box on
if nargin==5
saveas(gca,name)
end
end

function plotFigS2(abd1,Y1,expect1,Sc1,Sc2) % effect of non-functional group
cgabd=(groupsummary(abd1',expect1','sum'))';
[b,bint,~]=regress(Y1,x2fx(abd1));
Nstrains=size(abd1,2);

clf
W=33; H=36;
tiledlayout(2,2,"TileSpacing","compact")
set(gcf,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperSize',[W H],'PaperPosition',[0 0 W H], ...
     'Units','Centimeters','Position',[2 2 W H]);
ax1=nexttile;
plot(cgabd(:,2),Y1,'.')
xlabel("Abundance of non-functional group")
ylabel("Function")
h3=lsline;
h3.Color='r';
r=corr(Y1,cgabd(:,2));
text(2.6,0.4,"r = "+r,"FontSize",20)
set(gca,"FontSize",20)
box on
axis square

ax2=nexttile;
hold on
nfc=nnz(expect1==1);
% errorbar(0,b(1),b(1)-bint(1,1),bint(1,2)-b(1),'ko')
errorbar(1:nfc,b(2:(nfc+1)),b(2:(nfc+1))-bint(2:(nfc+1),1),bint(2:(nfc+1),2)-b(2:(nfc+1)),'ro')
errorbar(nfc+1:Nstrains,b((nfc+2):end),b((nfc+2):end)-bint((nfc+2):end,1),bint((nfc+2):end,2)-b((nfc+2):end),'bo')
hold off
yline(0,'--k')
% xticks(0:Nstrains)
xlabel("Species")
xlim([0 Nstrains+1])
ylabel("Coefficients")
set(gca,"FontSize",20)
box on
axis square

ax3=nexttile;
Sc1=mean(Sc1,3);
hold on
boxchart(Sc1,'MarkerSize',10,'MarkerColor',[0.3,0.3,0.3],'MarkerStyle','.','JitterOutliers','on','BoxEdgeColor',"k","BoxFaceColor",[0.5 0.5 0.5])
plot(1:5,mean(Sc1),'pentagramk','MarkerSize',12)
hold off
xticklabels(["EQO","EQO-2g","K-Means","Metropolis","Random"])
ylabel("Grouping score")
title("High competition")
ylim([0.2 1])
set(gca,"FontSize",20)
box on
axis square

ax4=nexttile;
Sc2=mean(Sc2,3);
hold on
boxchart(Sc2,'MarkerSize',10,'MarkerColor',[0.3,0.3,0.3],'MarkerStyle','.','JitterOutliers','on','BoxEdgeColor',"k","BoxFaceColor",[0.5 0.5 0.5])
plot(1:5,mean(Sc2),'pentagramk','MarkerSize',12)
hold off
xticklabels(["EQO","EQO-2g","K-Means","Metropolis","Random"])
ylabel("Grouping score")
title("Low competition")
ylim([0.2 1])
set(gca,"FontSize",20)
box on
axis square

pos= tightPosition(ax1,IncludeLabels=true);
x1=[pos(1)-0.01,pos(2)+pos(4),0.01,0.01]; x1=max(min(x1,1),0);
pos= tightPosition(ax2,IncludeLabels=true);
x2=[pos(1)-0.01,pos(2)+pos(4),0.01,0.01]; x2=max(min(x2,1),0);
pos= tightPosition(ax3,IncludeLabels=true);
x3=[pos(1)-0.01,pos(2)+pos(4)-0.01,0.01,0.01];x3=max(min(x3,1),0);
pos= tightPosition(ax4,IncludeLabels=true); 
x4=[pos(1)-0.01,pos(2)+pos(4)-0.01,0.01,0.01]; x4=max(min(x4,1),0);
annotation("textbox",x1,'String',"A","LineStyle",'none','FontSize',30)
annotation("textbox",x2,'String',"B","LineStyle",'none','FontSize',30)
annotation("textbox",x3,'String',"C","LineStyle",'none','FontSize',30)
annotation("textbox",x4,'String',"D","LineStyle",'none','FontSize',30)

exportgraphics(gcf,'figS2.pdf',"ContentType","vector")
end

function plotFigS4(GpE,GpK,GpM) % groupings for fig 2.
clf
W=50; H=40;
tiledlayout(2,3,"TileSpacing","tight","Padding","compact")
set(gcf,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperSize',[W H],'PaperPosition',[0 0 W H], ...
     'Units','Centimeters','Position',[2 2 W H]);
ax1=nexttile;
imagesc(GpE)
xlabel('Species')
ylabel('Datasets index')
title('EQO')
axis square
set(gca,"FontSize",25)

ax2=nexttile;
imagesc(GpK(:,:,1))
xlabel('Species')
ylabel('Datasets index')
title('2-group K-Means')
axis square
set(gca,"FontSize",25)

ax3=nexttile;
imagesc(GpM(:,:,1))
xlabel('Species')
ylabel('Datasets index')
title('2-group Metropolis')
axis square
set(gca,"FontSize",25)

ax4=nexttile(5);
imagesc(GpK(:,:,2))
xlabel('Species')
ylabel('Datasets index')
title('3-group K-Means')
axis square
set(gca,"FontSize",25)

ax5=nexttile(6);
imagesc(GpM(:,:,2))
xlabel('Species')
ylabel('Datasets index')
title('3-group Metropolis')
axis square
set(gca,"FontSize",25)

exportgraphics(gcf,'figS4.pdf',"ContentType","vector")
end

function plotFigS5(ScoreE,ScoreK,ScoreM,noise,Nsamples) % Heatmap of scores
me=mean(ScoreM,4);
km=mean(ScoreK,4);
eqo=mean(ScoreE,4);
clf
W=33; H=36;
tiledlayout(2,2,"TileSpacing","tight","Padding","compact")
set(gcf,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperSize',[W H],'PaperPosition',[0 0 W H], ...
     'Units','Centimeters','Position',[2 2 W H]);
ax1=nexttile;
imagesc(eqo)
colorbar
axis square xy
title('EQO')
set(gca,"FontSize",20)
yticks(1:length(Nsamples))
yticklabels(Nsamples)
ylabel('Number of samples')
xticks(1:length(noise))
xticklabels(noise)
xlabel('Relative noise')

ax2=nexttile;
imagesc(me(:,:,1))
colorbar
axis square xy
title('2-group Metropolis')
set(gca,"FontSize",20)
yticks(1:length(Nsamples))
yticklabels(Nsamples)
ylabel('Number of samples')
xticks(1:length(noise))
xticklabels(noise)
xlabel('Relative noise')

ax3=nexttile;
imagesc(km(:,:,1))
colorbar
axis square xy
title('2-group K-Means')
set(gca,"FontSize",20)
yticks(1:length(Nsamples))
yticklabels(Nsamples)
ylabel('Number of samples')
xticks(1:length(noise))
xticklabels(noise)
xlabel('Relative noise')

ax4=nexttile;
imagesc(km(:,:,2))
colorbar
axis square xy
title('3-group K-Means')
set(gca,"FontSize",20)
yticks(1:length(Nsamples))
yticklabels(Nsamples)
ylabel('Number of samples')
xticks(1:length(noise))
xticklabels(noise)
xlabel('Relative noise')

pos= tightPosition(ax1,IncludeLabels=true);
x1=[pos(1)-0.01,pos(2)+pos(4)-0.01,0.01,0.01]; x1=max(min(x1,1),0);
pos= tightPosition(ax2,IncludeLabels=true);
x2=[pos(1)-0.01,pos(2)+pos(4)-0.01,0.01,0.01]; x2=max(min(x2,1),0);
pos= tightPosition(ax3,IncludeLabels=true);
x3=[pos(1)-0.01,pos(2)+pos(4)-0.01,0.01,0.01]; x3=max(min(x3,1),0);
pos= tightPosition(ax4,IncludeLabels=true);
x4=[pos(1)-0.01,pos(2)+pos(4)-0.01,0.01,0.01]; x4=max(min(x4,1),0);
annotation("textbox",x1,'String',"A","LineStyle",'none','FontSize',30)
annotation("textbox",x2,'String',"B","LineStyle",'none','FontSize',30)
annotation("textbox",x3,'String',"C","LineStyle",'none','FontSize',30)
annotation("textbox",x4,'String',"D","LineStyle",'none','FontSize',30)

exportgraphics(gcf,'figS5.pdf',"ContentType","vector")
end