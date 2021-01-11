%-- 11/4/13, 10:23 AM --%
Mulfactau=0.33;
delta_h=18;
delta_h=-18;
v2=-110:0.5:40;
tauh=Mulfactau*(100/(7*exp((v2+60+delta_h)/11)+10*exp(-(v2+60+delta_h)/15))+0.6);
tauh=Mulfactau.*(100./(7.*exp((v2+60+delta_h)./11)+10.*exp(-(v2+60+delta_h)./15))+0.6);
plot(v2,tauh)
tauh2=Mulfactau.*(100./(7.*exp((v2+60)./11)+10.*exp(-(v2+60)./25))+0.6);
plot(v2,tauh2,'r')
hold on
plot(v2,tauh)
figure
plot(v2,tauh)
hold on
plot(v2,tauh2,'r')
tauh2=Mulfactau.*(100./(7.*exp((v2+60+18)./11)+10.*exp(-(v2+60+18)./25))+0.6);
plot(v2,tauh2,'r-')
tauh2=Mulfactau.*(100./(7.*exp((v2+60-18)./11)+10.*exp(-(v2+60-18)./25))+0.6);
plot(v2,tauh2,'r-')
%-- 11/4/13, 5:07 PM --%
cmp
doc cmp
test
silentsynapses_Gui
doc avona
PlotSilentSynapse
ss='/Users/lindameng/Documents/study/meng/102711/bs0003/SilentSynapseDistribute.eps';
pp='/Volumes/disk1/study'
doc regexprep
regexprep(ss,'/\wstudy',pp)
regexprep(ss,'/\w*y',pp)
regexprep(ss,'/\w*y$',pp)
regexprep(ss,'\w*y$',pp)
regexprep(ss,'\w*y\>',pp)
regexprep(ss,'/Users/lindameng/Documents/study',pp)
PlotSilentSynapse
avghistLayers
Cellrecord
Cellrecord(1)
Cellrecord{1}
guide
test
MeanValueCalculation
avghistLayers
guide
test
handles
test
handles.hd
hd.avgDensity
test
MeanValueCalculationSi
PK
Pk
find(~isnan(Pk))
hd.avgDensity
handles.hd
handles.Pchglayer
hd
handles
Zflag
find(Zflag>0)
Pk
find(Pk>0)
hd
hd.avgChgsi
Cellsum.avgChgsi
Cellsum.avgPk
handles.minmaps
doc scatterplot
doc statistic
doc scatte
guide
kmean
doc kmean
test
doc axis
test
doc plot
clear
Cellsum
plot(avgAgesi,avgPksi,'o')
plot(AvgAgesi,avgPksi,'o')
plot(Cellsum.AvgAgesi,Cellsum.avgPksi,'o')
plot(Cellsum.AvgAgesi,Cellsum.Pchglayer(:,2),'o')
plot(Cellsum.AvgAgesi,Cellsum.Pchglayer(:,1),'o')
plot(Cellsum.AvgAgesi,Cellsum.Pchglayer(:,3),'o')
hist(Cellsum.ratio)
doc hist
hist(Cellsum.ratio,10)
hist(Cellsum.ratio,20)
plot(Cellsum.AvgAgesi,Cellsum.ratio,'o')
doc hist
x=0:0.5:3;
hist(Cellsum.ratio,x)
x=0:0.1:3;
hist(Cellsum.ratio,x)
doc accum
cum
ccum
doc sum
doc cumsum
ecdf(Cellsum.ratio)
hold on
id=find(Cellsum.AvgAgesi>=3&Cellsum.AvgAgesi<6)
R35=Cellsum.ratio(id);
ecdf(R35)
id=find(Cellsum.AvgAgesi>=6&Cellsum.AvgAgesi<10)
R69=Cellsum.ratio(id);
ecdf(R69)
id=find(Cellsum.AvgAgesi>=10&Cellsum.AvgAgesi<16)
R10=Cellsum.ratio(id);
ecdf(R10)
hist(Cellsum.Pchglayer)
hist(Cellsum.Pchglayer[:,2])
hist(Cellsum.Pchglayer(:,2))
hist(Cellsum.Pchglayer(:,2),0:0.02:0.4)
hist(Cellsum.Pchglayer(:,2),0:0.02:0.5)
hist(Cellsum.Pchglayer(:,2),0:0.01:0.5)
figure(2)
hist(Cellsum.Ppeaklayer(:,2),0:0.01:0.5)
doc kmeans
[Id,c]=kmeans(Cellsum.Ppeaklayer(:,2), 2)
x=[Cellsum.Ppeaklayer(:,2),Cellsum.Pchgayer(:,2)];
x=[Cellsum.Ppeaklayer(:,2),Cellsum.Pchglayer(:,2)];
[Id,c]=kmeans(x, 2)
scatterhist(x,y,'Group',Id,'Location','SouthEast',...
'Direction','out','Color','kbr','LineStyle',{'-','-.',':'},...
'LineWidth',[2,2,2],'Marker','+od','MarkerSize',[4,5,6]);
scatterhist(x,'Group',Id,'Location','SouthEast',...
'Direction','out','Color','kbr','LineStyle',{'-','-.',':'},...
'LineWidth',[2,2,2],'Marker','+od','MarkerSize',[4,5,6]);
scatterhist(x(:,1),x(:,2),'Group',Id,'Location','SouthEast',...
'Direction','out','Color','kbr','LineStyle',{'-','-.',':'},...
'LineWidth',[2,2,2],'Marker','+od','MarkerSize',[4,5,6]);
figure(2)
hist(Cellsum.Ppeaklayer(:,2),0:0.01:0.5)
figure(3)
scatterhist(x(:,1),x(:,2),'Group',Id,'Location','SouthEast',...
'Direction','out','Color','kbr','LineStyle',{'-','-.',':'},...
'LineWidth',[2,2,2],'Marker','+od','MarkerSize',[4,5,6]);
[Id,c]=kmeans(x, 3)
scatterhist(x(:,1),x(:,2),'Group',Id,'Location','SouthEast',...
'Direction','out','Color','kbr','LineStyle',{'-','-.',':'},...
'LineWidth',[2,2,2],'Marker','+od','MarkerSize',[4,5,6]);
doc kmeans
D=[];
for j=1:100
[Id,C]=kmean(x,2);
D=[D; C(1,:)];
end
D2=[];
for j=1:100
[Id,C]=kmeans(x,2);
D=[D; C(1,:)];D2=[D2;C(2,:)];
end
x
idnan=find(~isnan(x(:,1)));
xx=x(idnan,:)
D=[];D2=[];for j=1:100
[Id,C]=kmean(xx,2);
D=[D; C(1,:)];D2=[D2; C(1,:)];
end
for j=1:100
[Id,C]=kmeans(xx,2);
D=[D; C(1,:)];D2=[D2; C(1,:)];
end
[Id,C]=kmeans(xx,2);
for j=1:100
[Id,C]=kmeans(xx,2);
D=[D; C(1,:)];D2=[D2; C(1,:)];
end
doc kmeans
[idx,ctrs] = kmeans(xx,2,...
'Distance','city',...
'Replicates',5,...
'Options',opts);
opts = statset('Display','final');
[idx,ctrs] = kmeans(xx,2,...
'Distance','city',...
'Replicates',5,...
'Options',opts);
for i=5:100
opts = statset('Display','final');
[idx,ctrs,sumd] = kmeans(xx,2,...
'Distance','city',...
'Replicates',i,...
'Options',opts);
D=[D,sumd];
emd
end
D=[];for i=5:100
[idx,ctrs,sumd] = kmeans(xx,2,...
'Distance','city',...
'Replicates',i,...
'Options',opts);
D=[D,sumd];
end
variance(D)
var(D)
D
size(D)
sumd
doc kmean
D=[];for i=5:100
[idx,ctrs,sumd] = kmeans(xx,2,...
'Distance','city',...
'Replicates',i);
D=[D,sumd];
end
var(D(1,:))
mean(D(1,:))
var(D(1,:))/mean(D(1,:))
[Id,C]=kmeans(xx,4);
[Id,C]=kmeans(xx,5);
[Id,C]=kmeans(xx,3);
[Id,C]=kmeans(xx,10);
[Id,C]=kmeans(xx,8);
[Id,C]=kmeans(xx,6);
[Id,C]=kmeans(xx,2);
figure(3)
scatterhist(xx(:,1),xx(:,2),'Group',Id,'Location','SouthEast',...
'Direction','out','Color','kbr','LineStyle',{'-','-.',':'},...
'LineWidth',[2,2,2],'Marker','+od','MarkerSize',[4,5,6]);
hold on
plot(C,'ko')
plot(C(1,:),'ko')
plot(C(:,1),'ko')
C
plot(C(1,:),'ko')
plot(C(1,1),C(1,2),'ko')
plot(C(2,1),C(2,2),'ko')
doc regression
plotregression(xx(:,1),xx(:,2))
%In plotregression at 107
D
Dchg=Cellsum.DistWidthChargeTp2HP60Drexp0Age0toInf;
x=Dchg(:,1);
y=Dchg(;,2)
y=Dchg(:,2);
X=[x y];
Id=kmeans(X,2);
scatterhist(x,y,'Group',Id,'Location','SouthEast',...
'Direction','out','Color','kbr','LineStyle',{'-','-.',':'},...
'LineWidth',[2,2,2],'Marker','+od','MarkerSize',[4,5,6]);
Dchg2=Cellsum.DistWidthChargeTp2HP50Drexp0Age0toInf;
x2=Dchg2(:,1);
y2=Dchg2(:,2);
scatterhist(x,y,'Group',Id,'Location','SouthEast',...
'Direction','out','Color','kbr','LineStyle',{'-','-.',':'},...
'LineWidth',[2,2,2],'Marker','+od','MarkerSize',[4,5,6]);
figure(2)
scatterhist(x,y,'Group',Id,'Location','SouthEast',...
'Direction','out','Color','kbr','LineStyle',{'-','-.',':'},...
'LineWidth',[2,2,2],'Marker','+od','MarkerSize',[4,5,6]);
scatterhist(x2,y2,'Group',Id,'Location','SouthEast',...
'Direction','out','Color','kbr','LineStyle',{'-','-.',':'},...
'LineWidth',[2,2,2],'Marker','+od','MarkerSize',[4,5,6]);
pca
doc pca
ddd
a(i).name
a
MeanValueCalculationSi
Cellrecord{1}
plot(Cellrecord{1}.L1chgTp2Hp60)
plot(Cellrecord{1}.L1pkTp2Hp60)
plot(Cellrecord{1}.L1PkTp2Hp60)
plot(Cellrecord{1}.L1PkTp2Hp60,'ko')
plot(Cellrecord{1}.L1chgTp2Hp60,'ko')
y1=Cellrecord{1}.L1chgTp2Hp60;
fn=1; Soma= Cellrecord{fn}.SomaCoordinates;
flipimg=Cellrecord{fn}.flipimg;
pth=char(Cellrecord{fn}.Pth);
stimcoordinates=Cellrecord{fn}.StimCoordinates;
rotate_angle=Cellrecord{fn}.SpatialRotation;
[stimCoordinates]=Matrixrotate_singlemap_GUI(Soma,stimcoordinates,flipimg,rotate_angle);
X=stimCoordinates(1,:);
Y=stimCoordinates(2,:);dx =Cellrecord{fn}.PatternSpacing(1); dy = dx;
xlim = [ceil(min(X))-dx-0.1 ceil(max(X))]; ylim = [ceil(min(Y))-dy-0.1 ceil(max(Y))];
xtick = xlim(1):dx:xlim(2); ytick = ylim(1):dy:ylim(2);
xi = xtick(1:end-1)+dx/2; yi = ytick(1:end-1)+dy/2;[mx ii]=max(y1);
jj=find(y1>mx*0.1);
if ~isempty(jj)
D1=abs(yi(jj(1))-yi(jj(end)));
else D1=0;
end
D1
y1
plot(yi,y1)
jj=find(y1>mx*0.1);
jj
mx
plot(yi,mx*0.1)
hold on
plot(yi,y1)
1*NaN
ScatterPlot_widthDistribution
axis image
ScatterPlot_widthDistribution
Cellrecord{fn}
dX=[dX,Cellrecord{fn}.PatternSpacing(1)];
Zflag=Cellrecord{fn}.meanflagPtxInd60;
IndZ=find(Zflag>0);
Zpk(IndZ)=Cellrecord{fn}.meanpeakPtxInd60(IndZ);dx =Cellrecord{fn}.PatternSpacing(1); dy = dx;
xlim = [ceil(min(X))-dx-0.1 ceil(max(X))]; ylim = [ceil(min(Y))-dy-0.1 ceil(max(Y))];
xtick = xlim(1):dx:xlim(2); ytick = ylim(1):dy:ylim(2);
xi = xtick(1:end-1)+dx/2; yi = ytick(1:end-1)+dy/2;
avg = gridavg2(X,Y,Zpk,xtick,ytick,false,1,0);
xx=nansum(avg);
yy=nansum(avg,2);
[mxx ii]=max(xx);
jj=find(xx>mxx*0.1);
if ~isempty(jj)
Dx=length(jj)*dx;
else Dx=0;
end
Zflag=Cellrecord{fn}.meanflagPtxInd60;
IndZ=find(Zflag>0);
Zpk(IndZ)=Cellrecord{fn}.meanpeakPtxInd60(IndZ);dx =Cellrecord{fn}.PatternSpacing(1); dy = dx;
xlim = [ceil(min(X))-dx-0.1 ceil(max(X))]; ylim = [ceil(min(Y))-dy-0.1 ceil(max(Y))];
xtick = xlim(1):dx:xlim(2); ytick = ylim(1):dy:ylim(2);
xi = xtick(1:end-1)+dx/2; yi = ytick(1:end-1)+dy/2;
avg = gridavg2(X,Y,Zpk,xtick,ytick,false,1,0);
xx=nansum(avg);
yy=nansum(avg,2);
[mxx ii]=max(xx);
jj=find(xx>mxx*0.1);
if ~isempty(jj)
Dx=length(jj)*dx;
else Dx=0;
end
Dx
[myy ii]=max(yy);
jj=find(yy>myy*0.1);
if ~isempty(jj)
Dy=length(jj)*dx;
else Dy=0;
end
Dy
ScatterPlot_widthDistribution
doc cumsum
ScatterPlot_widthDistribution
xi
xi(jj(end))
Bnd(4)
Bnd
xi(jj(1))-Bnd(4)
figure(100);plot(xi)
ScatterPlot_widthDistribution
doc ecdf
ScatterPlot_widthDistribution
doc ecdf
doc cdf
ScatterPlot_widthDistribution
size(X)
Y
squareform(Y)
T=squareform(Y)
image(T)
Z
image(Z)
c=cophenet(Z,Y)
Y = pdist(X,'cityblock');
Z = linkage(Y,'average');
c = cophenet(Z,Y)
cluster(Z)
cluster(Z,'cutoff',1.2)
Z
doc cluster
I=inconsistent(Z)
cluster(Z,'cutoff',2)
cluster(Z,'cutoff',3)
cluster(Z,'cutoff',4)
cluster(Z,'cutoff',5)
cluster(Z,'cutoff',10)
cluster(Z,'maxclust',2)
Id=cluster(Z,'maxclust',2)
scatterhist(y,x,'Group',Id,'Location','SouthEast',...
'Direction','out','Color','kbr','LineStyle',{'-','-.',':'},...
'LineWidth',[2,2,2],'Marker','+od','MarkerSize',[4,5,6]);
scatterhist(y,x,'Group',Id,'Location','SouthEast',...
'Direction','out','Color','kbr','LineStyle',{'-','-.',':'},...
'LineWidth',[2,2,2],'Marker','+od','MarkerSize',[4,5,6]);
scatterhist(y,xext,'Group',Id,'Location','SouthEast',...
'Direction','out','Color','kbr','LineStyle',{'-','-.',':'},...
'LineWidth',[2,2,2],'Marker','+od','MarkerSize',[4,5,6]);
Id=cluster(Z,'maxclust',3)
scatterhist(y,x,'Group',Id,'Location','SouthEast',...
'Direction','out','Color','kbr','LineStyle',{'-','-.',':'},...
'LineWidth',[2,2,2],'Marker','+od','MarkerSize',[4,5,6]);
scatterhist(y,xext,'Group',Id,'Location','SouthEast',...
'Direction','out','Color','kbr','LineStyle',{'-','-.',':'},...
'LineWidth',[2,2,2],'Marker','+od','MarkerSize',[4,5,6]);
Id=cluster(Z,'maxclust',4)
scatterhist(y,xext,'Group',Id,'Location','SouthEast',...
'Direction','out','Color','kbr','LineStyle',{'-','-.',':'},...
'LineWidth',[2,2,2],'Marker','+od','MarkerSize',[4,5,6]);
ScatterPlot_widthDistribution
doc marker
doc scatterhist
ScatterPlot_widthDistribution
PlotSilentSynapse
Cellrecord{j}
Cellrecord{1}
PlotSilentSynapse
axis ij
view 180
view (180,90)
view (270,90)
view (90,90)
PlotSilentSynapse
scatterdist
ScatterPlot_widthDistribution
figure(1)
close all
ScatterPlot_widthDistribution
figure(1)
ScatterPlot_widthDistribution
I=find(Id==1)
[a,b]=cedf(yy(I))
[a,b]=ecdf(yy(I))
figure(200)
plot(a,b)
[a,b]=ecdf(xext(I))
plot(a,b)
figure(200)
find(b>=0.89&&b<=0.9)
find(b>=0.89&b<=0.9)
find(b>=0.85&b<=0.9)
a(22)
plot(yy(I),xext(I),'.')
I=find(Id==2)
[a,b]=ecdf(xext(I))
plot(yy(I),xext(I),'.')
plot(a,b)
find(b>=0.85&b<=0.9)
find(b>=0.9&b<=0.92)
find(a>=0.9&a<=0.92)
find(a>=0.9&a<=0.91)
b(63)
find(a>=0.95&a<=0.951)
find(a>=0.95&a<=0.96)
b(67)
ScatterPlot_widthDistribution
mainClusterValidationNC(X)
X
size(X)
save('Dist.txt',X)
save('Dist.txt','X')
mainCVAP
X
doc fprintf
fprintf('Dist.txt','%f','X')
fid=fopen('Dist.txt','w');
fprintf(fid,'%f',X);
fprintf(fid,'%f/n',X);
fid=fopen('Dist.txt','w');
fprintf(fid,'%f\n',X);
X
fprintf(fid,'%f\n',X);
fid=fopen('Dist.txt','w');
fprintf(fid,'%f %f\n',X);
doc fprintf
fid=fopen('Dist.txt','w');
fprintf(fid,'%f\t%f\n',X');
R=[0 0.73205     0.88725     0.93041     0.96269      0.9711      0.9768     0.98252     0.98762     0.98972]
N=[1 2 3 4 5 6 7 8 9 10]
plot(N,R,'*-')
figure(100)
plot(N,R,'*-')
size(AvgAge)
size((Bnd(:,2)-Bnd(:,4))./(Bnd(:,1)-Bnd(:,4)))
X=[AvgAge',(Bnd(:,2)-Bnd(:,4))./(Bnd(:,1)-Bnd(:,4))];
ScatterPlot_widthDistribution
I=find(Id==2)
[a,b]=ecdf(xext(I));
find(a>=0.9&a<=0.91)
a(31)
b(31)
plot(a,b)
I=find(Id==1)
[a,b]=ecdf(xext(I));
find(a>=0.9&a<=0.91)
a(63)
b(63)
find(a>=0.91&a<=0.92)
b(64)
PlotSilentSynapse
view(-90,90)
view(90,0)
view(90,2700)
view(90,270)
view(90,180)
90
view(90,90)
doc vies
doc view
view(90,90,180)
doc view
doc rotate
rotate(gca,zdir,180)
rotate(gca,[1 0 0],180)
rotate(gca,[1 0 0],90)
rotate(gca,[1 1 0],90)
ScatterPlot_widthDistribution
Xaxis{1}
size(Xaxis{1})
ScatterPlot_widthDistribution
close all
ScatterPlot_widthDistribution
close all
ScatterPlot_widthDistribution
doc interplate
doc interp2
ScatterPlot_widthDistribution
[yi,xi]=meshgrid(yi,xi)
avg2=interp2(yi,xi,avg,yii,xii);
yi
size(yi)
size(xi)
size(avg)
[yii, xii]=meshgrid(yii,xii);
avg2=interp2(yi,xi,avg,yii,xii);
image(yii,xii,avg2)
doc interp2
doc surf
ScatterPlot_widthDistribution
doc
doc gridinterp
avg2=griddedInterpolant(yi,xi,avg,yii,xii)
doc griddedInterpolant
avg2=griddedInterpolant(yi,xi,avg)
size(yi)
size(xi)
size(avg)
doc
doc griddedInterpolant
[X,Y] = ndgrid(1:10,1:10);
X
[yi,xi]=meshgrid(yi,xi)
avg2=griddedInterpolant(yi,xi,avg)
doc griddedInterpolant
ScatterPlot_widthDistribution
[yi,xi]=ndgrid(yi,xi)
avg2=griddedInterpolant(yi,xi,avg)
xi
size(yi)
size(xi)
size(avg)
size(yi)
avg2=griddedInterpolant(yi,xi,avg')
doc griddedInterpolant
image(avg2)
imagesc(avg2)
avg2
imagesc(avg2.Values)
figure(100)
imagesc(avg2.Values)
size(avg2)
ScatterPlot_widthDistribution
[yi,xi]=meshgrid(yi,xi)
[yii, xii]=meshgrid(yii,xii);
avg2=interp2(yi,xi,avg,yii,xii);
imagesc(avg2)
doc interp2
avg2=interp2(yi,xi,avg,yii,xii,'cubic');
imagesc(avg2)
figure(100)
imagesc(avg2)
axis ij
imagesc(yi(1,:),xi(:,1),avg2)
avg2=interp2(yi,xi,avg,yii,xii);
imagesc(yi(1,:),xi(:,1),avg2)
ScatterPlot_widthDistribution
%-- 11/19/13, 7:55 PM --%
doc colormap
doc CLim
ScatterPlot_widthDistribution
doc plot
doc box
doc Axes
ScatterPlot_widthDistribution
close all
ScatterPlot_widthDistribution
doc ecdf
doc pareto
ScatterPlot_widthDistribution
doc
doc pareto
doc hist
doc bar
doc pareto
ScatterPlot_widthDistribution
bar(x,y)
doc cusum
yy=cumsum(y);
plot(x,yy)
hold on
bar(x,y)
ScatterPlot_widthDistribution
[Pchg60 Pchg50 Pchgsi]
ScatterPlot_widthDistribution
[Pchg60 Pchg50 Pchgsi]
ScatterPlot_widthDistribution
doc pam
ScatterPlot_widthDistribution
Indsi=find(Pchg60==0&Pchg50~=0);
hist(Pchg50(Indsi))
hist(Pchg50(Indsi),0:0.05:0.8)
hist(Pchg50(Indsi),0:0.02:0.8)
ScatterPlot_widthDistribution
Indsi=find(Pchg60==0&Pchg50~=0);
hist(Pchg50(Indsi),0:0.02:0.8)
hist(Pchg50(Indsi),0:0.02:0.4)
hist(Pchg50(Indsi),0:0.01:0.4)
Id1=Indsi
savepath= '/Volumes/disk1/silent';for i=1:length(Id1)
Xavg=[Xavg Xaxis{Id1(i)}];
Yavg=[Yavg Yaxis{Id1(i)}];
Zpk60=[Zpk60 Zpk60avg{Id1(i)}];
Zchg60=[Zchg60 Zchg60avg{Id1(i)}];
Zden60=[Zden60 Zden60avg{Id1(i)}];
Zpk50=[Zpk50 Zpk50avg{Id1(i)}];
Zchg50=[Zchg50 Zchg50avg{Id1(i)}];
Zden50=[Zden50 Zden50avg{Id1(i)}];
Zpksi=[Zpksi Zpksiavg{Id1(i)}];
Zchgsi=[Zchgsi Zchgsiavg{Id1(i)}];
Zdensi=[Zdensi Zdensiavg{Id1(i)}];
end
density=0.15;
dens=10;
N=length(Id1);
Ncell=N;
xlim = [-1000-dx-0.1 400];
ylim = [-650-dy-0.1 500];
xtick = xlim(1):dx:xlim(2);
ytick = ylim(1):dy:ylim(2);
xi = xtick(1:end-1)+dx/2;
yi = ytick(1:end-1)+dy/2;
xii=xi(1):1:xi(end);
yii=yi(1):0.5:yi(end);
% get rid of the places with few events
[cnt_exc_dir ind] =griddensity2(Xavg,Yavg,[],[],xtick,ytick,N,dens);
sel = (cnt_exc_dir<=dens);
isel = find(sel);
selrm = [];
for i=1:length(isel)
sel2 = find(ind==isel(i)); selrm = [selrm sel2];
end
Xavg(selrm) = [];
Zpk60(selrm) = [];
Zpksi(selrm) = [];
Zchg50(selrm) =[];
Zden60(selrm) = [];
Zdensi(selrm) = [];
Yavg(selrm) = [];
Zpk50(selrm)= [];
Zchg60(selrm) = [];
Zchgsi(selrm) = [];
Zden50(selrm) = [];
avg = gridavg2(Xavg,Yavg,Zden60,xtick,ytick,false,Ncell,0);
Ind2=find(avg<=density); avg(Ind2)=0;
currfig= find_figure(sprintf('Avgdensity60,Id==%i',c));
clf;
avg=avg';
[yi1,xi1]=meshgrid(yi,xi);
[yii,xii]=meshgrid(yii,xii);
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colormap(winter)
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
saveas(currfig,fullfile(savepath,sprintf('Avgdensity60Id%i.fig',c)),'fig');
avg = gridavg2(Xavg,Yavg,Zpk60,xtick,ytick,false,Ncell,0);
avg(Ind2)=0;
currfig= find_figure(sprintf('AvgPk60,Id==%i',c));
clf; avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
saveas(currfig,fullfile(savepath,sprintf('Avgpeak60Id%i.fig',c)),'fig');
avg = gridavg2(Xavg,Yavg,Zchg60,xtick,ytick,false,Ncell,0);
avg(Ind2)=0;
currfig= find_figure(sprintf('Avgchg60,Id==%i',c));
clf; avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
saveas(currfig,fullfile(savepath,sprintf('Avgchg60Id%i.fig',c)),'fig');
avg =gridavg2(Xavg,Yavg,Zden50,xtick,ytick,false,Ncell,0);
Ind50=find(avg<=density);
avg(Ind50)=0;
currfig= find_figure(sprintf('Avgdensity50_2,Id==%i',c));
clf;
avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colormap(winter)
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
saveas(currfig,fullfile(savepath,sprintf('Avgdensity50Id%i.fig',c)),'fig');
avg = gridavg2(Xavg,Yavg,Zpk50,xtick,ytick,false,Ncell,0); avg(Ind50)=0;
currfig= find_figure(sprintf('AvgPk50,Id==%i',c));
clf;
avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.51560693641618 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
saveas(currfig,fullfile(savepath,sprintf('Avgpeak50Id%i.fig',c)),'fig');
avg =gridavg2(Xavg,Yavg,Zchg50,xtick,ytick,false,Ncell,0);
avg(Ind50)=0;
currfig= find_figure(sprintf('Avgchg50,Id==%i',c));
clf;
avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 =axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
saveas(currfig,fullfile(savepath,sprintf('Avgchg50Id%i.fig',c)),'fig');
avg = gridavg2(Xavg,Yavg,Zdensi,xtick,ytick,false,Ncell,0);
Indsi=find(avg<=density); avg(Indsi)=0;
currfig= find_figure(sprintf('Avgdensitysi,Id==%i',c));
clf;
avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colormap(winter)
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
saveas(currfig,fullfile(savepath,sprintf('AvgdensitysiId%i.fig',c)),'fig');
avg= gridavg2(Xavg,Yavg,Zpksi,xtick,ytick,false,Ncell,0);
avg(Indsi)=0;
currfig= find_figure(sprintf('AvgPksi,Id==%i',c));
clf; avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
saveas(currfig,fullfile(savepath,sprintf('AvgpeaksiId%i.fig',c)),'fig');
avg =gridavg2(Xavg,Yavg,Zchgsi,xtick,ytick,false,Ncell,0);
avg(Indsi)=0;
currfig= find_figure(sprintf('Avgchgsi,Id==%i',c));
clf;
avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 =axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3=axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
saveas(currfig,fullfile(savepath,sprintf('AvgchgsiId%i.fig',c)),'fig');
c=1;
> savepath= '/Volumes/disk1/silent';for i=1:length(Id1)
Xavg=[Xavg Xaxis{Id1(i)}];
Yavg=[Yavg Yaxis{Id1(i)}];
Zpk60=[Zpk60 Zpk60avg{Id1(i)}];
Zchg60=[Zchg60 Zchg60avg{Id1(i)}];
Zden60=[Zden60 Zden60avg{Id1(i)}];
Zpk50=[Zpk50 Zpk50avg{Id1(i)}];
Zchg50=[Zchg50 Zchg50avg{Id1(i)}];
Zden50=[Zden50 Zden50avg{Id1(i)}];
Zpksi=[Zpksi Zpksiavg{Id1(i)}];
Zchgsi=[Zchgsi Zchgsiavg{Id1(i)}];
Zdensi=[Zdensi Zdensiavg{Id1(i)}];
end
density=0.15;
dens=10;
N=length(Id1);
Ncell=N;
xlim = [-1000-dx-0.1 400];
ylim = [-650-dy-0.1 500];
xtick = xlim(1):dx:xlim(2);
ytick = ylim(1):dy:ylim(2);
xi = xtick(1:end-1)+dx/2;
yi = ytick(1:end-1)+dy/2;
xii=xi(1):1:xi(end);
yii=yi(1):0.5:yi(end);
% get rid of the places with few events
[cnt_exc_dir ind] =griddensity2(Xavg,Yavg,[],[],xtick,ytick,N,dens);
sel = (cnt_exc_dir<=dens);
isel = find(sel);
selrm = [];
for i=1:length(isel)
sel2 = find(ind==isel(i)); selrm = [selrm sel2];
end
Xavg(selrm) = [];
Zpk60(selrm) = [];
Zpksi(selrm) = [];
Zchg50(selrm) =[];
Zden60(selrm) = [];
Zdensi(selrm) = [];
Yavg(selrm) = [];
Zpk50(selrm)= [];
Zchg60(selrm) = [];
Zchgsi(selrm) = [];
Zden50(selrm) = [];
avg = gridavg2(Xavg,Yavg,Zden60,xtick,ytick,false,Ncell,0);
Ind2=find(avg<=density); avg(Ind2)=0;
currfig= find_figure(sprintf('Avgdensity60,Id==%i',c));
clf;
avg=avg';
[yi1,xi1]=meshgrid(yi,xi);
[yii,xii]=meshgrid(yii,xii);
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colormap(winter)
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
saveas(currfig,fullfile(savepath,sprintf('Avgdensity60Id%i.fig',c)),'fig');
avg = gridavg2(Xavg,Yavg,Zpk60,xtick,ytick,false,Ncell,0);
avg(Ind2)=0;
currfig= find_figure(sprintf('AvgPk60,Id==%i',c));
clf; avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
saveas(currfig,fullfile(savepath,sprintf('Avgpeak60Id%i.fig',c)),'fig');
avg = gridavg2(Xavg,Yavg,Zchg60,xtick,ytick,false,Ncell,0);
avg(Ind2)=0;
currfig= find_figure(sprintf('Avgchg60,Id==%i',c));
clf; avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
saveas(currfig,fullfile(savepath,sprintf('Avgchg60Id%i.fig',c)),'fig');
avg =gridavg2(Xavg,Yavg,Zden50,xtick,ytick,false,Ncell,0);
Ind50=find(avg<=density);
avg(Ind50)=0;
currfig= find_figure(sprintf('Avgdensity50_2,Id==%i',c));
clf;
avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colormap(winter)
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
saveas(currfig,fullfile(savepath,sprintf('Avgdensity50Id%i.fig',c)),'fig');
avg = gridavg2(Xavg,Yavg,Zpk50,xtick,ytick,false,Ncell,0); avg(Ind50)=0;
currfig= find_figure(sprintf('AvgPk50,Id==%i',c));
clf;
avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.51560693641618 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
saveas(currfig,fullfile(savepath,sprintf('Avgpeak50Id%i.fig',c)),'fig');
avg =gridavg2(Xavg,Yavg,Zchg50,xtick,ytick,false,Ncell,0);
avg(Ind50)=0;
currfig= find_figure(sprintf('Avgchg50,Id==%i',c));
clf;
avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 =axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
saveas(currfig,fullfile(savepath,sprintf('Avgchg50Id%i.fig',c)),'fig');
avg = gridavg2(Xavg,Yavg,Zdensi,xtick,ytick,false,Ncell,0);
Indsi=find(avg<=density); avg(Indsi)=0;
currfig= find_figure(sprintf('Avgdensitysi,Id==%i',c));
clf;
avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colormap(winter)
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
saveas(currfig,fullfile(savepath,sprintf('AvgdensitysiId%i.fig',c)),'fig');
avg= gridavg2(Xavg,Yavg,Zpksi,xtick,ytick,false,Ncell,0);
avg(Indsi)=0;
currfig= find_figure(sprintf('AvgPksi,Id==%i',c));
clf; avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
saveas(currfig,fullfile(savepath,sprintf('AvgpeaksiId%i.fig',c)),'fig');
avg =gridavg2(Xavg,Yavg,Zchgsi,xtick,ytick,false,Ncell,0);
avg(Indsi)=0;
currfig= find_figure(sprintf('Avgchgsi,Id==%i',c));
clf;
avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 =axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3=axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
saveas(currfig,fullfile(savepath,sprintf('AvgchgsiId%i.fig',c)),'fig');
savepath= '/Volumes/disk1/silent';for i=1:length(Id1)
Xavg=[Xavg Xaxis{Id1(i)}];
Yavg=[Yavg Yaxis{Id1(i)}];
Zpk60=[Zpk60 Zpk60avg{Id1(i)}];
Zchg60=[Zchg60 Zchg60avg{Id1(i)}];
Zden60=[Zden60 Zden60avg{Id1(i)}];
Zpk50=[Zpk50 Zpk50avg{Id1(i)}];
Zchg50=[Zchg50 Zchg50avg{Id1(i)}];
Zden50=[Zden50 Zden50avg{Id1(i)}];
Zpksi=[Zpksi Zpksiavg{Id1(i)}];
Zchgsi=[Zchgsi Zchgsiavg{Id1(i)}];
Zdensi=[Zdensi Zdensiavg{Id1(i)}];
end
density=0.15;
dens=10;
N=length(Id1);
Ncell=N;
xlim = [-1000-dx-0.1 400];
ylim = [-650-dy-0.1 500];
xtick = xlim(1):dx:xlim(2);
ytick = ylim(1):dy:ylim(2);
xi = xtick(1:end-1)+dx/2;
yi = ytick(1:end-1)+dy/2;
xii=xi(1):1:xi(end);
yii=yi(1):0.5:yi(end);
% get rid of the places with few events
[cnt_exc_dir ind] =griddensity2(Xavg,Yavg,[],[],xtick,ytick,N,dens);
sel = (cnt_exc_dir<=dens);
isel = find(sel);
selrm = [];
for i=1:length(isel)
sel2 = find(ind==isel(i)); selrm = [selrm sel2];
end
Xavg(selrm) = [];
Zpk60(selrm) = [];
Zpksi(selrm) = [];
Zchg50(selrm) =[];
Zden60(selrm) = [];
Zdensi(selrm) = [];
Yavg(selrm) = [];
Zpk50(selrm)= [];
Zchg60(selrm) = [];
Zchgsi(selrm) = [];
Zden50(selrm) = [];
avg = gridavg2(Xavg,Yavg,Zden60,xtick,ytick,false,Ncell,0);
Ind2=find(avg<=density); avg(Ind2)=0;
currfig= find_figure(sprintf('Avgdensity60,Id==%i',c));
clf;
avg=avg';
[yi1,xi1]=meshgrid(yi,xi);
[yii,xii]=meshgrid(yii,xii);
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colormap(winter)
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
saveas(currfig,fullfile(savepath,sprintf('Avgdensity60Id%i.fig',c)),'fig');
avg = gridavg2(Xavg,Yavg,Zpk60,xtick,ytick,false,Ncell,0);
avg(Ind2)=0;
currfig= find_figure(sprintf('AvgPk60,Id==%i',c));
clf; avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
saveas(currfig,fullfile(savepath,sprintf('Avgpeak60Id%i.fig',c)),'fig');
avg = gridavg2(Xavg,Yavg,Zchg60,xtick,ytick,false,Ncell,0);
avg(Ind2)=0;
currfig= find_figure(sprintf('Avgchg60,Id==%i',c));
clf; avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
saveas(currfig,fullfile(savepath,sprintf('Avgchg60Id%i.fig',c)),'fig');
avg =gridavg2(Xavg,Yavg,Zden50,xtick,ytick,false,Ncell,0);
Ind50=find(avg<=density);
avg(Ind50)=0;
currfig= find_figure(sprintf('Avgdensity50_2,Id==%i',c));
clf;
avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colormap(winter)
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
saveas(currfig,fullfile(savepath,sprintf('Avgdensity50Id%i.fig',c)),'fig');
avg = gridavg2(Xavg,Yavg,Zpk50,xtick,ytick,false,Ncell,0); avg(Ind50)=0;
currfig= find_figure(sprintf('AvgPk50,Id==%i',c));
clf;
avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.51560693641618 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
saveas(currfig,fullfile(savepath,sprintf('Avgpeak50Id%i.fig',c)),'fig');
avg =gridavg2(Xavg,Yavg,Zchg50,xtick,ytick,false,Ncell,0);
avg(Ind50)=0;
currfig= find_figure(sprintf('Avgchg50,Id==%i',c));
clf;
avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 =axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
saveas(currfig,fullfile(savepath,sprintf('Avgchg50Id%i.fig',c)),'fig');
avg = gridavg2(Xavg,Yavg,Zdensi,xtick,ytick,false,Ncell,0);
Indsi=find(avg<=density); avg(Indsi)=0;
currfig= find_figure(sprintf('Avgdensitysi,Id==%i',c));
clf;
avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colormap(winter)
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
saveas(currfig,fullfile(savepath,sprintf('AvgdensitysiId%i.fig',c)),'fig');
avg= gridavg2(Xavg,Yavg,Zpksi,xtick,ytick,false,Ncell,0);
avg(Indsi)=0;
currfig= find_figure(sprintf('AvgPksi,Id==%i',c));
clf; avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
saveas(currfig,fullfile(savepath,sprintf('AvgpeaksiId%i.fig',c)),'fig');
avg =gridavg2(Xavg,Yavg,Zchgsi,xtick,ytick,false,Ncell,0);
avg(Indsi)=0;
currfig= find_figure(sprintf('Avgchgsi,Id==%i',c));
clf;
avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 =axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3=axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
saveas(currfig,fullfile(savepath,sprintf('AvgchgsiId%i.fig',c)),'fig');
Id1=find(Pchg60==0&Pchg50~=0);
Xavg=[];
Yavg=[];
dx=50;
dy=50;
Zpk60=[];
Zchg60=[];
Zden60=[];
Zpk50=[];
Zchg50=[];
Zden50=[];
Zpksi=[];
Zchgsi=[];
Zdensi=[];
for i=1:length(Id1)
Xavg=[Xavg Xaxis{Id1(i)}];
Yavg=[Yavg Yaxis{Id1(i)}];
Zpk60=[Zpk60 Zpk60avg{Id1(i)}];
Zchg60=[Zchg60 Zchg60avg{Id1(i)}];
Zden60=[Zden60 Zden60avg{Id1(i)}];
Zpk50=[Zpk50 Zpk50avg{Id1(i)}];
Zchg50=[Zchg50 Zchg50avg{Id1(i)}];
Zden50=[Zden50 Zden50avg{Id1(i)}];
Zpksi=[Zpksi Zpksiavg{Id1(i)}];
Zchgsi=[Zchgsi Zchgsiavg{Id1(i)}];
Zdensi=[Zdensi Zdensiavg{Id1(i)}];
end
density=0.15;
dens=10;
N=length(Id1);
Ncell=N;
xlim = [-1000-dx-0.1 400];
ylim = [-650-dy-0.1 500];
xtick = xlim(1):dx:xlim(2);
ytick = ylim(1):dy:ylim(2);
xi = xtick(1:end-1)+dx/2;
yi = ytick(1:end-1)+dy/2;
xii=xi(1):1:xi(end);
yii=yi(1):0.5:yi(end);
% get rid of the places with few events
[cnt_exc_dir ind] =griddensity2(Xavg,Yavg,[],[],xtick,ytick,N,dens);
sel = (cnt_exc_dir<=dens);
isel = find(sel);
selrm = [];
for i=1:length(isel)
sel2 = find(ind==isel(i)); selrm = [selrm sel2];
end
Xavg(selrm) = [];
Zpk60(selrm) = [];
Zpksi(selrm) = [];
Zchg50(selrm) =[];
Zden60(selrm) = [];
Zdensi(selrm) = [];
Yavg(selrm) = [];
Zpk50(selrm)= [];
Zchg60(selrm) = [];
Zchgsi(selrm) = [];
Zden50(selrm) = [];
avg = gridavg2(Xavg,Yavg,Zden60,xtick,ytick,false,Ncell,0);
Ind2=find(avg<=density); avg(Ind2)=0;
currfig= find_figure(sprintf('Avgdensity60,Id==%i',c));
clf;
avg=avg';
[yi1,xi1]=meshgrid(yi,xi);
[yii,xii]=meshgrid(yii,xii);
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colormap(winter)
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
avg = gridavg2(Xavg,Yavg,Zpk60,xtick,ytick,false,Ncell,0);
avg(Ind2)=0;
currfig= find_figure(sprintf('AvgPk60,Id==%i',c));
clf; avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
saveas(currfig,fullfile(savepath,sprintf('Avgpeak60Id%i.fig',c)),'fig');
avg = gridavg2(Xavg,Yavg,Zchg60,xtick,ytick,false,Ncell,0);
avg(Ind2)=0;
currfig= find_figure(sprintf('Avgchg60,Id==%i',c));
clf; avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colorbar
title(sprintf('#cell %i',Ncell))
hold on
currfig= find_figure(sprintf('Avgdensity50_2,Id==%i',c));
clf;
avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colormap(winter)
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
saveas(currfig,fullfile(savepath,sprintf('Avgdensity50Id%i.fig',c)),'fig');
avg = gridavg2(Xavg,Yavg,Zpk50,xtick,ytick,false,Ncell,0); avg(Ind50)=0;
currfig= find_figure(sprintf('AvgPk50,Id==%i',c));
clf;
avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colorbar
title(sprintf('#cell %i',Ncell))
hold on
avg =gridavg2(Xavg,Yavg,Zden50,xtick,ytick,false,Ncell,0);
Ind50=find(avg<=density);
avg(Ind50)=0;
currfig= find_figure(sprintf('Avgdensity50_2,Id==%i',c));
clf;
avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colormap(winter)
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
dens=0.1;
avg =gridavg2(Xavg,Yavg,Zden50,xtick,ytick,false,Ncell,0);
Ind50=find(avg<=density);
avg(Ind50)=0;
currfig= find_figure(sprintf('Avgdensity50_2,Id==%i',c));
clf;
avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colormap(winter)
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
density=0.1;
avg =gridavg2(Xavg,Yavg,Zden50,xtick,ytick,false,Ncell,0);
Ind50=find(avg<=density);
avg(Ind50)=0;
currfig= find_figure(sprintf('Avgdensity50_2,Id==%i',c));
clf;
avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colormap(winter)
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
avg = gridavg2(Xavg,Yavg,Zpk50,xtick,ytick,false,Ncell,0); avg(Ind50)=0;
currfig= find_figure(sprintf('AvgPk50,Id==%i',c));
clf;
avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.51560693641618 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
avg = gridavg2(Xavg,Yavg,Zchg60,xtick,ytick,false,Ncell,0);
avg(Ind2)=0;
currfig= find_figure(sprintf('Avgchg60,Id==%i',c));
clf; avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
xlim = [ceil(min(Xavg))-dx-0.1 ceil(max(Xavg))]; ylim = [ceil(min(Yavg))-dy-0.1 ceil(max(Yavg))];
xtick = xlim(1):dx:xlim(2); ytick = ylim(1):dy:ylim(2);
xi = xtick(1:end-1)+dx/2; yi = ytick(1:end-1)+dy/2;
avg = gridavg2(Xavg,Yavg,Zchg60,xtick,ytick,false,Ncell,0);
avg(Ind2)=0;
currfig= find_figure(sprintf('Avgchg60,Id==%i',c));
clf; avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
xii=xi(1):1:xi(end);
yii=yi(1):0.5:yi(end);
[yi1,xi1]=meshgrid(yi,xi);
[yii,xii]=meshgrid(yii,xii);
avg = gridavg2(Xavg,Yavg,Zchg60,xtick,ytick,false,Ncell,0);
avg(Ind2)=0;
currfig= find_figure(sprintf('Avgchg60,Id==%i',c));
clf; avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
avg = gridavg2(Xavg,Yavg,Zden60,xtick,ytick,false,Ncell,0);
Ind2=find(avg<=density); avg(Ind2)=0;
currfig= find_figure(sprintf('Avgdensity60,Id==%i',c));
clf;
avg=avg';
[yi1,xi1]=meshgrid(yi,xi);
[yii,xii]=meshgrid(yii,xii);
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colormap(winter)
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
%-- 11/20/13, 7:40 PM --%
ScatterPlot_widthDistribution
xlim = [ceil(min(Xavg))-dx-0.1 ceil(max(Xavg))]; ylim = [ceil(min(Yavg))-dy-0.1 ceil(max(Yavg))];
xtick = xlim(1):dx:xlim(2); ytick = ylim(1):dy:ylim(2);
xi = xtick(1:end-1)+dx/2; yi = ytick(1:end-1)+dy/2;
Id1=find(Pchg60==0&Pchg50~=0);
Xavg=[];
Yavg=[];
dx=50;
dy=50;
Zpk60=[];
Zchg60=[];
Zden60=[];
Zpk50=[];
Zchg50=[];
Zden50=[];
Zpksi=[];
Zchgsi=[];
Zdensi=[];
for i=1:length(Id1)
Xavg=[Xavg Xaxis{Id1(i)}];
Yavg=[Yavg Yaxis{Id1(i)}];
Zpk60=[Zpk60 Zpk60avg{Id1(i)}];
Zchg60=[Zchg60 Zchg60avg{Id1(i)}];
Zden60=[Zden60 Zden60avg{Id1(i)}];
Zpk50=[Zpk50 Zpk50avg{Id1(i)}];
Zchg50=[Zchg50 Zchg50avg{Id1(i)}];
Zden50=[Zden50 Zden50avg{Id1(i)}];
Zpksi=[Zpksi Zpksiavg{Id1(i)}];
Zchgsi=[Zchgsi Zchgsiavg{Id1(i)}];
Zdensi=[Zdensi Zdensiavg{Id1(i)}];
end
xlim = [ceil(min(Xavg))-dx-0.1 ceil(max(Xavg))]; ylim = [ceil(min(Yavg))-dy-0.1 ceil(max(Yavg))];
xtick = xlim(1):dx:xlim(2); ytick = ylim(1):dy:ylim(2);
xi = xtick(1:end-1)+dx/2; yi = ytick(1:end-1)+dy/2;
xii=xi(1):1:xi(end);
yii=yi(1):0.5:yi(end);
[yi1,xi1]=meshgrid(yi,xi);
[yii,xii]=meshgrid(yii,xii);
currfig= find_figure(sprintf('Avgdensity60,Id==%i',c));
clf;
avg=avg';
[yi1,xi1]=meshgrid(yi,xi);
[yii,xii]=meshgrid(yii,xii);
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colormap(winter)
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
c=1;
currfig= find_figure(sprintf('Avgdensity60,Id==%i',c));
clf;
avg=avg';
[yi1,xi1]=meshgrid(yi,xi);
[yii,xii]=meshgrid(yii,xii);
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colormap(winter)
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
avg = gridavg2(Xavg,Yavg,Zden60,xtick,ytick,false,Ncell,0);
Ind2=find(avg<=density); avg(Ind2)=0;
currfig= find_figure(sprintf('Avgdensity60,Id==%i',c));
clf;
avg=avg';
[yi1,xi1]=meshgrid(yi,xi);
[yii,xii]=meshgrid(yii,xii);
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colormap(winter)
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
Ncell=N
Ncell=Nsi
avg = gridavg2(Xavg,Yavg,Zden60,xtick,ytick,false,Ncell,0);
Ind2=find(avg<=density); avg(Ind2)=0;
currfig= find_figure(sprintf('Avgdensity60,Id==%i',c));
clf;
avg=avg';
[yi1,xi1]=meshgrid(yi,xi);
[yii,xii]=meshgrid(yii,xii);
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colormap(winter)
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
density=0.1;
avg = gridavg2(Xavg,Yavg,Zden60,xtick,ytick,false,Ncell,0);
Ind2=find(avg<=density); avg(Ind2)=0;
currfig= find_figure(sprintf('Avgdensity60,Id==%i',c));
clf;
avg=avg';
[yi1,xi1]=meshgrid(yi,xi);
[yii,xii]=meshgrid(yii,xii);
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colormap(winter)
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
%-- 11/20/13, 7:52 PM --%
ScatterPlot_widthDistribution
xlim = [ceil(min(Xavg))-dx-0.1 ceil(max(Xavg))]; ylim = [ceil(min(Yavg))-dy-0.1 ceil(max(Yavg))];
xtick = xlim(1):dx:xlim(2); ytick = ylim(1):dy:ylim(2);
xi = xtick(1:end-1)+dx/2; yi = ytick(1:end-1)+dy/2;
Id1=find(Pchg60==0&Pchg50~=0);
Xavg=[];
Yavg=[];
dx=50;
dy=50;
Zpk60=[];
Zchg60=[];
Zden60=[];
Zpk50=[];
Zchg50=[];
Zden50=[];
Zpksi=[];
Zchgsi=[];
Zdensi=[];
for i=1:length(Id1)
Xavg=[Xavg Xaxis{Id1(i)}];
Yavg=[Yavg Yaxis{Id1(i)}];
Zpk60=[Zpk60 Zpk60avg{Id1(i)}];
Zchg60=[Zchg60 Zchg60avg{Id1(i)}];
Zden60=[Zden60 Zden60avg{Id1(i)}];
Zpk50=[Zpk50 Zpk50avg{Id1(i)}];
Zchg50=[Zchg50 Zchg50avg{Id1(i)}];
Zden50=[Zden50 Zden50avg{Id1(i)}];
Zpksi=[Zpksi Zpksiavg{Id1(i)}];
Zchgsi=[Zchgsi Zchgsiavg{Id1(i)}];
Zdensi=[Zdensi Zdensiavg{Id1(i)}];
end
Id1=find(Pchg60==0&Pchg50~=0);
Xavg=[];
Yavg=[];
dx=50;
dy=50;
Zpk60=[];
Zchg60=[];
Zden60=[];
Zpk50=[];
Zchg50=[];
Zden50=[];
Zpksi=[];
Zchgsi=[];
Zdensi=[];
for i=1:length(Id1)
Xavg=[Xavg Xaxis{Id1(i)}];
Yavg=[Yavg Yaxis{Id1(i)}];
Zpk60=[Zpk60 Zpk60avg{Id1(i)}];
Zchg60=[Zchg60 Zchg60avg{Id1(i)}];
Zden60=[Zden60 Zden60avg{Id1(i)}];
Zpk50=[Zpk50 Zpk50avg{Id1(i)}];
Zchg50=[Zchg50 Zchg50avg{Id1(i)}];
Zden50=[Zden50 Zden50avg{Id1(i)}];
Zpksi=[Zpksi Zpksiavg{Id1(i)}];
Zchgsi=[Zchgsi Zchgsiavg{Id1(i)}];
Zdensi=[Zdensi Zdensiavg{Id1(i)}];
end
density=0.15;
dens=10;
N=length(Id1);
Ncell=N;
xii=xi(1):1:xi(end);
yii=yi(1):0.5:yi(end);
xlim = [ceil(min(Xavg))-dx-0.1 ceil(max(Xavg))]; ylim = [ceil(min(Yavg))-dy-0.1 ceil(max(Yavg))];
xtick = xlim(1):dx:xlim(2); ytick = ylim(1):dy:ylim(2);
xi = xtick(1:end-1)+dx/2; yi = ytick(1:end-1)+dy/2;
xi = xtick(1:end-1)+dx/2;
yi = ytick(1:end-1)+dy/2;
xii=xi(1):1:xi(end);
yii=yi(1):0.5:yi(end);
[cnt_exc_dir ind] =griddensity2(Xavg,Yavg,[],[],xtick,ytick,N,dens);
sel = (cnt_exc_dir<=dens);
isel = find(sel);
selrm = [];
for i=1:length(isel)
sel2 = find(ind==isel(i)); selrm = [selrm sel2];
end
Xavg(selrm) = [];
Zpk60(selrm) = [];
Zpksi(selrm) = [];
Zchg50(selrm) =[];
Zden60(selrm) = [];
Zdensi(selrm) = [];
Yavg(selrm) = [];
Zpk50(selrm)= [];
Zchg60(selrm) = [];
Zchgsi(selrm) = [];
Zden50(selrm) = [];
[yi1,xi1]=meshgrid(yi,xi);
[yii,xii]=meshgrid(yii,xii);
avg = gridavg2(Xavg,Yavg,Zden60,xtick,ytick,false,Ncell,0);
Ind2=find(avg<=density); avg(Ind2)=0;
currfig= find_figure(sprintf('Avgdensity60,Id==%i',c));
clf;
avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colormap(winter)
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
c=1
[cnt_exc_dir ind] =griddensity2(Xavg,Yavg,[],[],xtick,ytick,N,dens);
sel = (cnt_exc_dir<=dens);
isel = find(sel);
selrm = [];
for i=1:length(isel)
sel2 = find(ind==isel(i)); selrm = [selrm sel2];
end
Xavg(selrm) = [];
Zpk60(selrm) = [];
Zpksi(selrm) = [];
Zchg50(selrm) =[];
Zden60(selrm) = [];
Zdensi(selrm) = [];
Yavg(selrm) = [];
Zpk50(selrm)= [];
Zchg60(selrm) = [];
Zchgsi(selrm) = [];
Zden50(selrm) = [];
avg = gridavg2(Xavg,Yavg,Zden60,xtick,ytick,false,Ncell,0);
Ind2=find(avg<=density); avg(Ind2)=0;
currfig= find_figure(sprintf('Avgdensity60,Id==%i',c));
clf;
avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colormap(winter)
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
avg =gridavg2(Xavg,Yavg,Zden50,xtick,ytick,false,Ncell,0);
Ind50=find(avg<=density);
avg(Ind50)=0;
currfig= find_figure(sprintf('Avgdensity50_2,Id==%i',c));
clf;
avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colormap(winter)
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
density=0.1;
avg =gridavg2(Xavg,Yavg,Zden50,xtick,ytick,false,Ncell,0);
Ind50=find(avg<=density);
avg(Ind50)=0;
currfig= find_figure(sprintf('Avgdensity50_2,Id==%i',c));
clf;
avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colormap(winter)
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
Id1=find(Id==1);
Id2=find(Pchg60(Id1)==0&Pchg50(Id1)~=0);
Xavg=[];
Yavg=[];
dx=50;
dy=50;
Zpk60=[];
Zchg60=[];
Zden60=[];
Zpk50=[];
Zchg50=[];
Zden50=[];
Zpksi=[];
Zchgsi=[];
Zdensi=[];
for i=1:length(Id2)
Xavg=[Xavg Xaxis{Id1(Id2(i))}];
Yavg=[Yavg Yaxis{Id1(Id2(i))}];
Zpk60=[Zpk60 Zpk60avg{Id1(Id2(i))}];
Zchg60=[Zchg60 Zchg60avg{Id1(Id2(i))}];
Zden60=[Zden60 Zden60avg{Id1(Id2(i))}];
Zpk50=[Zpk50 Zpk50avg{Id1(Id2(i))}];
Zchg50=[Zchg50 Zchg50avg{Id1(Id2(i))}];
Zden50=[Zden50 Zden50avg{Id1(Id2(i))}];
Zpksi=[Zpksi Zpksiavg{Id1(Id2(i))}];
Zchgsi=[Zchgsi Zchgsiavg{Id1(Id2(i))}];
Zdensi=[Zdensi Zdensiavg{Id1(Id2(i))}];
end
[cnt_exc_dir ind] =griddensity2(Xavg,Yavg,[],[],xtick,ytick,N,dens);
sel = (cnt_exc_dir<=dens);
isel = find(sel);
selrm = [];
for i=1:length(isel)
sel2 = find(ind==isel(i)); selrm = [selrm sel2];
end
Xavg(selrm) = [];
Zpk60(selrm) = [];
Zpksi(selrm) = [];
Zchg50(selrm) =[];
Zden60(selrm) = [];
Zdensi(selrm) = [];
Yavg(selrm) = [];
Zpk50(selrm)= [];
Zchg60(selrm) = [];
Zchgsi(selrm) = [];
Zden50(selrm) = [];
Ind2=find(avg<=density); avg(Ind2)=0;
currfig= find_figure(sprintf('Avgdensity60,Id==%i',c));
clf;
avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colormap(winter)
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
[cnt_exc_dir ind] =griddensity2(Xavg,Yavg,[],[],xtick,ytick,N,dens);
sel = (cnt_exc_dir<=dens);
isel = find(sel);
selrm = [];
for i=1:length(isel)
sel2 = find(ind==isel(i)); selrm = [selrm sel2];
end
Xavg(selrm) = [];
Zpk60(selrm) = [];
Zpksi(selrm) = [];
Zchg50(selrm) =[];
Zden60(selrm) = [];
Zdensi(selrm) = [];
Yavg(selrm) = [];
Zpk50(selrm)= [];
Zchg60(selrm) = [];
Zchgsi(selrm) = [];
Zden50(selrm) = [];
size(yi1)
avg = gridavg2(Xavg,Yavg,Zden60,xtick,ytick,false,Ncell,0);
Ind2=find(avg<=density); avg(Ind2)=0;
currfig= find_figure(sprintf('Avgdensity60,Id==%i',c));
clf;
avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colormap(winter)
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
density=.1
avg = gridavg2(Xavg,Yavg,Zden60,xtick,ytick,false,Ncell,0);
Ind2=find(avg<=density); avg(Ind2)=0;
currfig= find_figure(sprintf('Avgdensity60,Id==%i',c));
clf;
avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colormap(winter)
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
Id1=find(Id==2);
Id2=find(Pchg60(Id1)==0&Pchg50(Id1)~=0);
Xavg=[];
Yavg=[];
dx=50;
dy=50;
Zpk60=[];
Zchg60=[];
Zden60=[];
Zpk50=[];
Zchg50=[];
Zden50=[];
Zpksi=[];
Zchgsi=[];
Zdensi=[];
for i=1:length(Id2)
Xavg=[Xavg Xaxis{Id1(Id2(i))}];
Yavg=[Yavg Yaxis{Id1(Id2(i))}];
Zpk60=[Zpk60 Zpk60avg{Id1(Id2(i))}];
Zchg60=[Zchg60 Zchg60avg{Id1(Id2(i))}];
Zden60=[Zden60 Zden60avg{Id1(Id2(i))}];
Zpk50=[Zpk50 Zpk50avg{Id1(Id2(i))}];
Zchg50=[Zchg50 Zchg50avg{Id1(Id2(i))}];
Zden50=[Zden50 Zden50avg{Id1(Id2(i))}];
Zpksi=[Zpksi Zpksiavg{Id1(Id2(i))}];
Zchgsi=[Zchgsi Zchgsiavg{Id1(Id2(i))}];
Zdensi=[Zdensi Zdensiavg{Id1(Id2(i))}];
end
xlim = [ceil(min(Xavg))-dx-0.1 ceil(max(Xavg))]; ylim = [ceil(min(Yavg))-dy-0.1 ceil(max(Yavg))];
xtick = xlim(1):dx:xlim(2); ytick = ylim(1):dy:ylim(2);
xi = xtick(1:end-1)+dx/2; yi = ytick(1:end-1)+dy/2;
xi = xtick(1:end-1)+dx/2;
yi = ytick(1:end-1)+dy/2;
xii=xi(1):1:xi(end);
yii=yi(1):0.5:yi(end);
[cnt_exc_dir ind] =griddensity2(Xavg,Yavg,[],[],xtick,ytick,N,dens);
sel = (cnt_exc_dir<=dens);
isel = find(sel);
selrm = [];
for i=1:length(isel)
sel2 = find(ind==isel(i)); selrm = [selrm sel2];
end
Xavg(selrm) = [];
Zpk60(selrm) = [];
Zpksi(selrm) = [];
Zchg50(selrm) =[];
Zden60(selrm) = [];
Zdensi(selrm) = [];
Yavg(selrm) = [];
Zpk50(selrm)= [];
Zchg60(selrm) = [];
Zchgsi(selrm) = [];
Zden50(selrm) = [];
Ind2=find(avg<=density); avg(Ind2)=0;
[yi1,xi1]=meshgrid(yi,xi);
[yii,xii]=meshgrid(yii,xii);
[cnt_exc_dir ind] =griddensity2(Xavg,Yavg,[],[],xtick,ytick,N,dens);
sel = (cnt_exc_dir<=dens);
isel = find(sel);
selrm = [];
for i=1:length(isel)
sel2 = find(ind==isel(i)); selrm = [selrm sel2];
end
Xavg(selrm) = [];
Zpk60(selrm) = [];
Zpksi(selrm) = [];
Zchg50(selrm) =[];
Zden60(selrm) = [];
Zdensi(selrm) = [];
Yavg(selrm) = [];
Zpk50(selrm)= [];
Zchg60(selrm) = [];
Zchgsi(selrm) = [];
Zden50(selrm) = [];
avg = gridavg2(Xavg,Yavg,Zden60,xtick,ytick,false,Ncell,0);
Ind2=find(avg<=density); avg(Ind2)=0;
currfig= find_figure(sprintf('Avgdensity60,Id==%i',c));
clf;
avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colormap(winter)
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
avg =gridavg2(Xavg,Yavg,Zden50,xtick,ytick,false,Ncell,0);
Ind50=find(avg<=density);
avg(Ind50)=0;
currfig= find_figure(sprintf('Avgdensity50_2,Id==%i',c));
clf;
avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colormap(winter)
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
avg =gridavg2(Xavg,Yavg,Zden50,xtick,ytick,false,Ncell,0);
Ind50=find(avg<=density);
avg(Ind50)=0;
currfig= find_figure(sprintf('Avgdensity50_2,Id==%i',c));
clf;
avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colormap(winter)
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
density=0.08
avg =gridavg2(Xavg,Yavg,Zden50,xtick,ytick,false,Ncell,0);
Ind50=find(avg<=density);
avg(Ind50)=0;
currfig= find_figure(sprintf('Avgdensity50_2,Id==%i',c));
clf;
avg=avg';
avg2=interp2(yi1,xi1,avg,yii,xii);
axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
'DataAspectRatio',[1 1 1]);
imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
colormap(winter)
colorbar
title(sprintf('#cell %i',Ncell))
hold on
plot(0,0,'wo','MarkerSize',15)
hold(axes1,'all');
xmat=sum(avg);
ymat=sum(avg,2);
axes2 = axes('Parent',currfig,'XAxisLocation','top',...
'Position',[0.228901734104046 0.00151975683889925 0.515606936416185 0.106382978723411]);
plot(yi,xmat,'k-','LineWidth',2)
set(gca,'XLim',[min(ytick) max(ytick)])
box off
axes3 = axes('Parent',currfig,'YDir','reverse',...
'Position',[0.134086260560249 0.109118541033435 0.0961538461538463 0.813373860182371]);
plot(xi,ymat,'k-','LineWidth',2)
set(gca,'XLim',[min(xtick) max(xtick)])
view(90,90)
ScatterPlot_widthDistribution
close all
ScatterPlot_widthDistribution
close all
ScatterPlot_widthDistribution
Id1=find(Id==1&Ag<=5&Ag>=3);
Id2=find(Id==1&Ag<=9&Ag>=6);
Id3=find(Id==1&Ag<=15&Ag>=10);
P1=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=9&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=10))]
Id1=find(Id==2&Ag<=5&Ag>=3);
Id2=find(Id==2&Ag<=9&Ag>=6);
Id3=find(Id==2&Ag<=15&Ag>=10);
P2=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=9&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=10))]
xx=[1 2 3]
plot(xx,P1,'r*')
figure(100)
plot(xx,P1,'r*')
plot(xx,P1,'r*-')
hold on
plot(xx,P2,'b*-')
9
Id1=find(Id==1&Ag<=6&Ag>=3);
Id2=find(Id==1&Ag<=9&Ag>=7);
Id3=find(Id==1&Ag<=15&Ag>=10);
P1=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=9&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=10))]
Id1=find(Id==2&Ag<=6&Ag>=3);
Id2=find(Id==2&Ag<=9&Ag>=7);
Id3=find(Id==2&Ag<=15&Ag>=10);
P1=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=9&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=10))]
P1=[length(Id1)/length(find(Ag<=6&Ag>=3)) length(Id2)/length(find(Ag<=9&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=10))]
Id1=find(Id==1&Ag<=5&Ag>=3);
Id2=find(Id==1&Ag<=10&Ag>=6);
Id3=find(Id==1&Ag<=15&Ag>=11);
P1=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=10&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=11))]
Id1=find(Id==1&Ag<=5&Ag>=3);
Id2=find(Id==1&Ag<=10&Ag>=6);
Id3=find(Id==1&Ag<=15&Ag>=11);
P1=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=10&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=11))]
Id1=find(Id==1&Ag<=5&Ag>=3);
Id2=find(Id==1&Ag<=10&Ag>=6);
Id3=find(Id==1&Ag<=15&Ag>=11);
P1=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=10&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=11))]
Id1=find(Id==1&Ag<=6&Ag>=3);
Id2=find(Id==1&Ag<=10&Ag>=7);
Id3=find(Id==1&Ag<=15&Ag>=11);
P1=[length(Id1)/length(find(Ag<=6&Ag>=3)) length(Id2)/length(find(Ag<=10&Ag>=7)) length(Id3)/length(find(Ag<=15&Ag>=11))]
Id1=find(Id==1&Ag<=5&Ag>=3);
Id2=find(Id==1&Ag<=8&Ag>=6);
Id3=find(Id==1&Ag<=15&Ag>=9);
P1=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=8&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=9))]
plot(xx,P1,'b*-')
plot(xx,1-P1,'r*-')
length(find(Id==1))
ScatterPlot_widthDistribution
Id1=find((Pk60>0.1||Pk50>0.1)&Ag<=5&Ag>=3);
Id2=find((Pk60>0.1||Pk50>0.1)&Ag<=8&Ag>=6);
Id3=find((Pk60>0.1||Pk50>0.1)&Ag<=15&Ag>=9);
P1=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=8&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=9))]
Id1=find((Pk60>0.1||Pk50>0.1)&Ag<=5&Ag>=3);
Id2=find((Ppk60>0.1||Ppk50>0.1)&Ag<=8&Ag>=6);
Id3=find((Ppk60>0.1||Ppk50>0.1)&Ag<=15&Ag>=9);
P1=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=8&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=9))]
Id1=find((Ppk60>0.1||Ppk50>0.1)&Ag<=5&Ag>=3);
Id2=find((Ppk60>0.1||Ppk50>0.1)&Ag<=8&Ag>=6);
Id3=find((Ppk60>0.1||Ppk50>0.1)&Ag<=15&Ag>=9);
P1=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=8&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=9))]
Id1=find((Ppk60>0.1|Ppk50>0.1)&Ag<=5&Ag>=3);
Id2=find((Ppk60>0.1|Ppk50>0.1)&Ag<=8&Ag>=6);
Id3=find((Ppk60>0.1|Ppk50>0.1)&Ag<=15&Ag>=9);
P1=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=8&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=9))]
Id1=find((Ppk60>0.1|Ppk50>0.1)&Ag<=5&Ag>=3);
Id2=find((Ppk60>0.1|Ppk50>0.1)&Ag<=9&Ag>=6);
Id3=find((Ppk60>0.1|Ppk50>0.1)&Ag<=15&Ag>=10);
P1=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=9&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=10))]
X=[Ppk60,yy]
Id=kmeans(X,clusterN);
scatterhist(xext2, Ppk60,'Group',Id,'Location','SouthEast',...
%     'Direction','out','Color','kbrm','LineStyle',{'-','-.',':','-'},...
%     'LineWidth',[2,2,2,2],'Marker','+od*','MarkerSize',[4,5,6,7]);
scatterhist(xext2, Ppk60,'Group',Id,'Location','SouthEast',...    'Direction','out','Color','kbrm','LineStyle',{'-','-.',':','-'},...
%     'LineWidth',[2,2,2,2],'Marker','+od*','MarkerSize',[4,5,6,7]);
scatterhist(xext, Ppk60,'Group',Id,'Location','SouthEast',...
'Direction','out','Color','kbrm','LineStyle',{'-','-.',':','-'},...
'LineWidth',[2,2,2,2],'Marker','+od*','MarkerSize',[4,5,6,7]);
X
X=[Ppk60,xext]
Id=kmeans(X,clusterN);
scatterhist(xext2, Ppk60,'Group',Id,'Location','SouthEast',...
'Direction','out','Color','kbrm','LineStyle',{'-','-.',':','-'},...
'LineWidth',[2,2,2,2],'Marker','+od*','MarkerSize',[4,5,6,7]);
scatterhist(xext, Ppk60,'Group',Id,'Location','SouthEast',...
'Direction','out','Color','kbrm','LineStyle',{'-','-.',':','-'},...
'LineWidth',[2,2,2,2],'Marker','+od*','MarkerSize',[4,5,6,7]);
X=[Ppk60,xext]
Id=kmeans(X,3);
scatterhist(xext2, Ppk60,'Group',Id,'Location','SouthEast',...
'Direction','out','Color','kbrm','LineStyle',{'-','-.',':','-'},...
'LineWidth',[2,2,2,2],'Marker','+od*','MarkerSize',[4,5,6,7]);
scatterhist(xext, Ppk60,'Group',Id,'Location','SouthEast',...
'Direction','out','Color','kbrm','LineStyle',{'-','-.',':','-'},...
'LineWidth',[2,2,2,2],'Marker','+od*','MarkerSize',[4,5,6,7]);
scatterhist(xext2, Ppk60+Ppeaklayer60(:,1),'Group',Id,'Location','SouthEast',...
'Direction','out','Color','kbrm','LineStyle',{'-','-.',':','-'},...
'LineWidth',[2,2,2,2],'Marker','+od*','MarkerSize',[4,5,6,7]);
scatterhist(xext, Ppk60+Ppeaklayer60(:,1),'Group',Id,'Location','SouthEast',...
'Direction','out','Color','kbrm','LineStyle',{'-','-.',':','-'},...
'LineWidth',[2,2,2,2],'Marker','+od*','MarkerSize',[4,5,6,7]);
Id=kmeans(X,2);
scatterhist(xext, Ppk60+Ppeaklayer60(:,1),'Group',Id,'Location','SouthEast',...
'Direction','out','Color','kbrm','LineStyle',{'-','-.',':','-'},...
'LineWidth',[2,2,2,2],'Marker','+od*','MarkerSize',[4,5,6,7]);
close all
scatterhist(xext, Ppk60+Ppeaklayer60(:,1),'Group',Id,'Location','SouthEast',...
'Direction','out','Color','kbrm','LineStyle',{'-','-.',':','-'},...
'LineWidth',[2,2,2,2],'Marker','+od*','MarkerSize',[4,5,6,7]);
ScatterPlot_widthDistribution
close all
ScatterPlot_widthDistribution
close all
ScatterPlot_widthDistribution
X=[yy Ppk60];
Id=kmeans(X,clusterN);
figure(2)
scatterhist(yy, Ppk60,'Group',Id,'Location','SouthEast',...
'Direction','out','Color','kbrm','LineStyle',{'-','-.',':','-'},...
'LineWidth',[2,2,2,2],'Marker','+od*','MarkerSize',[4,5,6,7]);
ScatterPlot_widthDistribution
close all
ScatterPlot_widthDistribution
test
guide                                                                                                   gu
guide
a
ScatterPlot_widthDistribution
folderMeanname
folderMeanpath
mapave_gui2
guide
%-- 11/22/13, 2:01 PM --%
test
guide
ScatterPlot_widthDistribution
close all
ScatterPlot_widthDistribution
fopen('Dist.txt','w')
h=fopen('Dist.txt','w');
doc fopen
fprintf(h,'%f\t%f\n',X')
fclose(h)
X
mainCVAP
ScatterPlot_widthDistribution
fopen('Dist.txt','w')
h=fopen('Dist.txt','w');
doc fopen
fprintf(h,'%f\t%f\n',X')
fclose(h)
X
ind=find(isnan(Ppk60));
Ppk60(ind)=0;
X=[xext Ppk60];
fopen('Dist1.txt','w')
h=fopen('Dist1.txt','w');
doc fopen
fprintf(h,'%f\t%f\n',X')
fclose(h)
R=[0 0.70832     0.91146     0.95618     0.98064     0.98064     0.98064     0.98571 ]
tt=0:1:10
figure(100) plot(tt, R)
figure(100); plot(tt, R)
R=[0 0.70832     0.91146     0.95618     0.98064     0.98064     0.98064     0.98571      0.9888     0.98794]
figure(100); plot(tt, R)
size(tt)
size(R)
tt=1:1:10
figure(100); plot(tt, R)
figure(100); plot(tt, R,'r*-','LineWidth','2')
figure(100); plot(tt, R,'r*-','LineWidth',2)
figure(100); plot(tt, R,'r*-','LineWidth',1)
Id1=find(Id==1&Ag<=5&Ag>=3);
Id2=find(Id==1&Ag<=8&Ag>=6);
Id3=find(Id==1&Ag<=15&Ag>=9);
P1=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=8&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=9))]
size(find(Id==1))
size(find(Id==2))
figure(100);plot([1,2,3],P1)
figure(100);plot([1,2,3],P1,'r*-','LineWidth','1')
figure(100);plot([1,2,3],P1,'r*-','LineWidth',1)
hold on
plot([1,2,3],1-P1,'k*-','LineWidth',1)
ScatterPlot_widthDistribution
close all
ScatterPlot_widthDistribution
hist(y(Id1))
figure(1)
hist(y(Id1),0:50:1000)
xhist=0:50:1000; figure(1)
subplot(2,3,1)
[a,b]=hist(y,xhist);
a=a./sum(a);
bar(b,a)
hold on
yy=cumsum(a);
plot(b,yy)
subplot(2,3,4)
[a,b]=hist(y(Id2),xhist)
a=a./sum(a);
bar(b,a)
hold on
yy=cumsum(a);
plot(b,yy)
subplot(2,3,2)
[a,b]=hist(y2(Id1),xhist)
a=a./sum(a);
bar(b,a)
hold on
yy=cumsum(a);
plot(b,yy)
subplot(2,3,5)
[a,b]=hist(y2(Id2),xhist)
y=y./sum(y);
bar(x,y)
hold on
yy=cumsum(y);
plot(x,yy)
subplot(2,3,3)
[a,b]=hist(y3(Id1),xhist)
a=a./sum(a);
bar(b,a)
hold on
yy=cumsum(a);
plot(b,yy)
subplot(2,3,6)
[y,x]=hist(y3(Id2),xhist)
a=a./sum(a);
bar(b,a)
hold on
yy=cumsum(a);
plot(b,yy)
xhist=0:50:1000; figure(1)
subplot(2,3,1)
[a,b]=hist(y,xhist);
a=a./sum(a);
bar(b,a)
hold on
yy=cumsum(a);
plot(b,yy)
subplot(2,3,4)
[a,b]=hist(y(Id2),xhist)
a=a./sum(a);
bar(b,a)
hold on
yy=cumsum(a);
plot(b,yy)
subplot(2,3,2)
[a,b]=hist(y2(Id1),xhist)
a=a./sum(a);
bar(b,a)
hold on
yy=cumsum(a);
plot(b,yy)
subplot(2,3,5)
[a,b]=hist(y2(Id2),xhist)
a=a./sum(a);
bar(b,a)
hold on
yy=cumsum(a);
plot(x,yy)
subplot(2,3,3)
[a,b]=hist(y3(Id1),xhist)
a=a./sum(a);
bar(b,a)
hold on
yy=cumsum(a);
plot(b,yy)
subplot(2,3,6)
[y,x]=hist(y3(Id2),xhist)
a=a./sum(a);
bar(b,a)
hold on
yy=cumsum(a);
plot(b,yy)
subplot(2,3,1)
[a,b]=hist(y,xhist);
a=a./sum(a);
bar(b,a)
hold on
yy=cumsum(a);
plot(b,yy)
subplot(2,3,4)
[a,b]=hist(y(Id2),xhist)
a=a./sum(a);
bar(b,a)
hold on
yy=cumsum(a);
plot(b,yy)
subplot(2,3,2)
[a,b]=hist(y2(Id1),xhist)
a=a./sum(a);
bar(b,a)
hold on
yy=cumsum(a);
plot(b,yy)
subplot(2,3,5)
[a,b]=hist(y2(Id2),xhist)
a=a./sum(a);
bar(b,a)
hold on
yy=cumsum(a);
plot(b,yy)
subplot(2,3,3)
[a,b]=hist(y3(Id1),xhist)
a=a./sum(a);
bar(b,a)
hold on
yy=cumsum(a);
plot(b,yy)
subplot(2,3,6)
[y,x]=hist(y3(Id2),xhist)
a=a./sum(a);
bar(b,a)
hold on
yy=cumsum(a);
plot(b,yy)
ecdf(y(Id1))
size(Id1)
doc ecdf
figure(2)
ecdf(y(Id1))
ScatterPlot_widthDistribution
ecdf(y(Id1))
xhist=0:40:1000;
ecdf(y(Id1),xhist)
doc ecdf
ecdf(y(Id1))
doc cumsum
Y=[];
for i=1:length(xhist)
t=y(Id1)>=xhist(i)
end
for i=2:length(xhist)
t=y(Id1)>=xhist(i-1)&y(Id1)<xhist(i);
Y=[Y sum(t)];
end
plot(xhist+20,Y)
size(xhist)
size(Y)
plot(xhist(1:end-1)+20,Y)
Y=[];
for i=2:length(xhist)
t=y(y(Id1)<=xhist(i);
Y=[Y sum(t)];
end
for i=2:length(xhist)
t=y(Id1)<=xhist(i);
Y=[Y sum(t)];
end
plot(xhist(1:end-1)+20,Y)
Y=[];
for i=2:length(xhist)
t=y(Id2)<=xhist(i);
Y=[Y sum(t)];
end
hold on
plot(xhist(1:end-1)+20,Y)
Y=[];
for i=2:length(xhist)
t=y(Id1)<=xhist(i);
Y=[Y sum(t)];
end
Y=Y./length(Id1);
plot(xhist(1:end-1)+20,Y)
Y=[];
for i=2:length(xhist)
t=y(Id2)<=xhist(i);
Y=[Y sum(t)];
end
Y=Y./length(Id2);
hold on
plot(xhist(1:end-1)+20,Y)
plot(xhist(1:end-1)+20,Y,'r')
Y=[];
for i=1:length(xhist)
t=y(Id2)<=xhist(i);
Y=[Y sum(t)];
end
Y=Y./length(Id2);
p
plot(xhist+20,Y,'r')
ScatterPlot_widthDistribution
Y=[];
for i=1:length(xhist)
t=y(Id2)<=xhist(i);
size(xhist)
end
size(xhist)
ScatterPlot_widthDistribution
mC=[nanmean(y(Id1)) nanmean(y2(Id1)) nanmean(y3(Id1))];
vC=[nanstd(y(Id1)) nanstd(y2(Id1)) nanstd(y3(Id1))];
mC=[mC;nanmean(y(Id2)) nanmean(y2(Id2)) nanmean(y3(Id2))];
vC=[vC;nanstd(y(Id2)) nanstd(y2(Id2)) nanstd(y3(Id2))]
bar(mC)
ScatterPlot_widthDistribution
hist(y(Id1),0:50:1000)
bar(mC)
size(Id1)
ScatterPlot_widthDistribution
bar(mC)
hold on
doc error
doc errorbar
xa=[0.8 1 1.2; 1.8 2 2.2];
errorbar(xa,mC,vC)
bar(mC)
hold on
xa=[0.75 1 1.25; 1.75 2 2.25];
errorbar(xa,mC,vC,'k*')
bar(mC)
hold on
xa=[0.78 1 1.22; 1.78 2 2.22];
errorbar(xa,mC,vC,'k*')
doc boxplot
C=[y(Id1) y(Id2)];
C=[y(Id1); y(Id2)];
G=[ones(length(Id1),1);ones(length(Id2),1)];
G=[ones(length(Id1),1);ones(length(Id2),1).*2];
boxplot(C,G)
ranksum(y(Id1),y(Id2))
ranksum(y2(Id1),y2(Id2))
ranksum(y3(Id1),y3(Id2))
ScatterPlot_widthDistribution
%-- 11/25/13, 10:04 AM --%
ScatterPlot_widthDistribution
close all
ScatterPlot_widthDistribution
close all
ScatterPlot_widthDistribution
close all
ScatterPlot_widthDistribution
close all
ScatterPlot_widthDistribution
Id1=find(Id==1&Ag<=5&Ag>=3);
Id2=find(Id==1&Ag<=9&Ag>=6);
Id3=find(Id==1&Ag<=15&Ag>=10);
P1=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=9&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=10))]
size(find(Id==1))
size(find(Id==2))
figure(100);plot([1,2,3],P1,'r*-','LineWidth',1)
hold on
plot([1,2,3],1-P1,'k*-','LineWidth',1)
Id1=find(Id==1&Ag<=5&Ag>=3);
Id2=find(Id==1&Ag<=8Ag>=6);
Id3=find(Id==1&Ag<=15&Ag>=9);
P1=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=8&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=9))]
size(find(Id==1))
size(find(Id==2))
figure(100);plot([1,2,3],P1,'r*-','LineWidth',1)
hold on
plot([1,2,3],1-P1,'k*-','LineWidth',1)
Id1=find(Id==1&Ag<=5&Ag>=3);
Id2=find(Id==1&Ag<=8&Ag>=6);
Id3=find(Id==1&Ag<=15&Ag>=9);
P1=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=8&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=9))]
size(find(Id==1))
size(find(Id==2))
figure(100);plot([1,2,3],P1,'r*-','LineWidth',1)
hold on
plot([1,2,3],1-P1,'k*-','LineWidth',1)
Id1=find(Ppk60>=0.1&Ag<=5&Ag>=3);
Id2=find(Ppk60>=0.1&Ag<=8Ag>=6);
Id3=find(Ppk60>=0.1&Ag<=15&Ag>=9);
P1=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=8&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=9))]
size(find(Id==1))
size(find(Id==2))
figure(100);plot([1,2,3],P1,'r*-','LineWidth',1)
hold on
plot([1,2,3],1-P1,'k*-','LineWidth',1)
Id1=find(Ppk60>=0.1&Ag<=5&Ag>=3);
Id2=find(Ppk60>=0.1&Ag<=8&Ag>=6);
Id3=find(Ppk60>=0.1&Ag<=15&Ag>=9);
P1=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=8&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=9))]
size(find(Id==1))
size(find(Id==2))
figure(100);plot([1,2,3],P1,'r*-','LineWidth',1)
hold on
plot([1,2,3],1-P1,'k*-','LineWidth',1)
Id1=find(Ppk50>=0.1&Ag<=5&Ag>=3);
Id2=find(Ppk50>=0.1&Ag<=8&Ag>=6);
Id3=find(Ppk50>=0.1&Ag<=15&Ag>=9);
P1=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=8&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=9))]
size(find(Id==1))
size(find(Id==2))
figure(100);plot([1,2,3],P1,'r*-','LineWidth',1)
hold on
plot([1,2,3],1-P1,'k*-','LineWidth',1)
Id1=find(Ppk50>=0.05&Ag<=5&Ag>=3);
Id2=find(Ppk50>=0.05&Ag<=8&Ag>=6);
Id3=find(Ppk50>=0.05&Ag<=15&Ag>=9);
P1=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=8&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=9))]
size(find(Id==1))
size(find(Id==2))
figure(100);plot([1,2,3],P1,'r*-','LineWidth',1)
hold on
plot([1,2,3],1-P1,'k*-','LineWidth',1)
close all;Id1=find(Ppk60>=0.01&Ag<=5&Ag>=3);
Id2=find(Ppk60>=0.01&Ag<=8&Ag>=6);
Id3=find(Ppk60>=0.01&Ag<=15&Ag>=9);
P1=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=8&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=9))]
size(find(Id==1))
size(find(Id==2))
figure(100);plot([1,2,3],P1,'r*-','LineWidth',1)
hold on
plot([1,2,3],1-P1,'k*-','LineWidth',1)
ScatterPlot_widthDistribution
close all;Id1=find(Ppk60>=0.01&Ag<=5&Ag>=3);
Id2=find(Ppk60>=0.01&Ag<=8&Ag>=6);
Id3=find(Ppk60>=0.01&Ag<=15&Ag>=9);
P1=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=8&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=9))]
size(find(Id==1))
size(find(Id==2))
figure(100);plot([1,2,3],P1,'r*-','LineWidth',1)
hold on
plot([1,2,3],1-P1,'k*-','LineWidth',1)
close all;Id1=find(Ppk60>=0.01&Ag<=5&Ag>=3);
Id2=find(Ppk60>=0.02&Ag<=8&Ag>=6);
Id3=find(Ppk60>=0.02&Ag<=15&Ag>=9);
P1=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=8&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=9))]
size(find(Id==1))
size(find(Id==2))
figure(100);plot([1,2,3],P1,'r*-','LineWidth',1)
hold on
plot([1,2,3],1-P1,'k*-','LineWidth',1)
close all;Id1=find(Ppk60>=0.03&Ag<=5&Ag>=3);
Id2=find(Ppk60>=0.03&Ag<=8&Ag>=6);
Id3=find(Ppk60>=0.03&Ag<=15&Ag>=9);
P1=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=8&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=9))]
size(find(Id==1))
size(find(Id==2))
figure(100);plot([1,2,3],P1,'r*-','LineWidth',1)
hold on
plot([1,2,3],1-P1,'k*-','LineWidth',1)
close all;Id1=find(Ppk60>=0.04&Ag<=5&Ag>=3);
Id2=find(Ppk60>=0.04&Ag<=8&Ag>=6);
Id3=find(Ppk60>=0.04&Ag<=15&Ag>=9);
P1=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=8&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=9))]
size(find(Id==1))
size(find(Id==2))
figure(100);plot([1,2,3],P1,'r*-','LineWidth',1)
hold on
plot([1,2,3],1-P1,'k*-','LineWidth',1)
close all;Id1=find(Ppk50>=0.01&Ag<=5&Ag>=3);
Id2=find(Ppk50>=0.01&Ag<=8&Ag>=6);
Id3=find(Ppk50>=0.01&Ag<=15&Ag>=9);
P1=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=8&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=9))]
size(find(Id==1))
size(find(Id==2))
figure(100);plot([1,2,3],P1,'r*-','LineWidth',1)
hold on
plot([1,2,3],1-P1,'k*-','LineWidth',1)
close all;Id1=find(Ppk50>=0.02&Ag<=5&Ag>=3);
Id2=find(Ppk50>=0.02&Ag<=8&Ag>=6);
Id3=find(Ppk50>=0.02&Ag<=15&Ag>=9);
P1=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=8&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=9))]
size(find(Id==1))
size(find(Id==2))
figure(100);plot([1,2,3],P1,'r*-','LineWidth',1)
hold on
plot([1,2,3],1-P1,'k*-','LineWidth',1)
close all;Id1=find(Ppk50>=0.02&Ag<=5&Ag>=3);
Id2=find(Ppk50>=0.02&Ag<=8&Ag>=6);
Id3=find(Ppk50>=0.02&Ag<=15&Ag>=9);
P1=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=8&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=9))]
size(find(Id==1))
size(find(Id==2))
figure(100);plot([1,2,3],P1,'r*-','LineWidth',1)
hold on
plot([1,2,3],1-P1,'k*-','LineWidth',1)
close all;Id1=find(Ppk50>=0.03&Ag<=5&Ag>=3);
Id2=find(Ppk50>=0.03&Ag<=8&Ag>=6);
Id3=find(Ppk50>=0.03&Ag<=15&Ag>=9);
P1=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=8&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=9))]
size(find(Id==1))
size(find(Id==2))
figure(100);plot([1,2,3],P1,'r*-','LineWidth',1)
hold on
plot([1,2,3],1-P1,'k*-','LineWidth',1)
close all;Id1=find(Ppk50>=0.01&Ag<=5&Ag>=3);
Id2=find(Ppk50>=0.01&Ag<=8&Ag>=6);
Id3=find(Ppk50>=0.01&Ag<=15&Ag>=9);
P1=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=8&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=9))]
size(find(Id==1))
size(find(Id==2))
figure(100);plot([1,2,3],P1,'r*-','LineWidth',1)
hold on
plot([1,2,3],1-P1,'k*-','LineWidth',1)
close all;Id1=find(Ppk50>=0.05&Ag<=5&Ag>=3);
Id2=find(Ppk50>=0.05&Ag<=8&Ag>=6);
Id3=find(Ppk50>=0.05&Ag<=15&Ag>=9);
P1=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=8&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=9))]
size(find(Id==1))
size(find(Id==2))
figure(100);plot([1,2,3],P1,'r*-','LineWidth',1)
hold on
plot([1,2,3],1-P1,'k*-','LineWidth',1)
close all;Id1=find(Ppk60>=0.05&Ag<=5&Ag>=3);
Id2=find(Ppk60>=0.05&Ag<=8&Ag>=6);
Id3=find(Ppk60>=0.05&Ag<=15&Ag>=9);
P1=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=8&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=9))]
size(find(Id==1))
size(find(Id==2))
figure(100);plot([1,2,3],P1,'r*-','LineWidth',1)
hold on
plot([1,2,3],1-P1,'k*-','LineWidth',1)
close all;Id1=find(Ppk60>=0.02&Ag<=5&Ag>=3);
Id2=find(Ppk60>=0.02&Ag<=8&Ag>=6);
Id3=find(Ppk60>=0.02&Ag<=15&Ag>=9);
P1=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=8&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=9))]
size(find(Id==1))
size(find(Id==2))
figure(100);plot([1,2,3],P1,'r*-','LineWidth',1)
hold on
plot([1,2,3],1-P1,'k*-','LineWidth',1)
close all;Id1=find(Ppk60>=0.02&Ag<=5&Ag>=3);
Id2=find(Ppk60>=0.02&Ag<=9&Ag>=6);
Id3=find(Ppk60>=0.02&Ag<=15&Ag>=10);
P1=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=9&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=10))]
size(find(Id==1))
size(find(Id==2))
figure(100);plot([1,2,3],P1,'r*-','LineWidth',1)
hold on
plot([1,2,3],1-P1,'k*-','LineWidth',1)
close all;Id1=find(Ppk60>=0.05&Ag<=5&Ag>=3);
Id2=find(Ppk60>=0.05&Ag<=9&Ag>=6);
Id3=find(Ppk60>=0.05&Ag<=15&Ag>=10);
P1=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=9&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=10))]
size(find(Id==1))
size(find(Id==2))
figure(100);plot([1,2,3],P1,'r*-','LineWidth',1)
hold on
plot([1,2,3],1-P1,'k*-','LineWidth',1)
close all;Id1=find(Ppk60>=0.04&Ag<=5&Ag>=3);
Id2=find(Ppk60>=0.04&Ag<=9&Ag>=6);
Id3=find(Ppk60>=0.04&Ag<=15&Ag>=10);
P1=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=9&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=10))]
size(find(Id==1))
size(find(Id==2))
figure(100);plot([1,2,3],P1,'r*-','LineWidth',1)
hold on
plot([1,2,3],1-P1,'k*-','LineWidth',1)
close all;Id1=find(Ppk60>=0.03&Ag<=5&Ag>=3);
Id2=find(Ppk60>=0.03&Ag<=9&Ag>=6);
Id3=find(Ppk60>=0.03&Ag<=15&Ag>=10);
P1=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=9&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=10))]
size(find(Id==1))
size(find(Id==2))
figure(100);plot([1,2,3],P1,'r*-','LineWidth',1)
hold on
plot([1,2,3],1-P1,'k*-','LineWidth',1)
close all;Id1=find(Ppk60>=0.02&Ag<=5&Ag>=3);
Id2=find(Ppk60>=0.02&Ag<=9&Ag>=6);
Id3=find(Ppk60>=0.02&Ag<=15&Ag>=10);
P1=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=9&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=10))]
size(find(Id==1))
size(find(Id==2))
figure(100);plot([1,2,3],P1,'r*-','LineWidth',1)
hold on
plot([1,2,3],1-P1,'k*-','LineWidth',1)
close all;Id1=find(Ppk60>=0.01&Ag<=5&Ag>=3);
Id2=find(Ppk60>=0.01&Ag<=9&Ag>=6);
Id3=find(Ppk60>=0.01&Ag<=15&Ag>=10);
P1=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=9&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=10))]
size(find(Id==1))
size(find(Id==2))
figure(100);plot([1,2,3],P1,'r*-','LineWidth',1)
hold on
plot([1,2,3],1-P1,'k*-','LineWidth',1)
close all;Id1=find(Ppk60>=0.005&Ag<=5&Ag>=3);
Id2=find(Ppk60>=0.005&Ag<=9&Ag>=6);
Id3=find(Ppk60>=0.005&Ag<=15&Ag>=10);
P1=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=9&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=10))]
size(find(Id==1))
size(find(Id==2))
figure(100);plot([1,2,3],P1,'r*-','LineWidth',1)
hold on
plot([1,2,3],1-P1,'k*-','LineWidth',1)
x=0:0.005:1;
hist(Ppk60,x)
x=0:0.05:1;
hist(Ppk60,x)
x=0:0.01:1;
hist(Ppk60,x)
close all;Id1=find(Ppk60>=0.02&Ag<=5&Ag>=3);
Id2=find(Ppk60>=0.02&Ag<=9&Ag>=6);
Id3=find(Ppk60>=0.02&Ag<=15&Ag>=10);
P1=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=9&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=10))]
size(find(Id==1))
size(find(Id==2))
figure(100);plot([1,2,3],P1,'r*-','LineWidth',1)
hold on
plot([1,2,3],1-P1,'k*-','LineWidth',1)
mean(Ppk60(Id1))
mean(Ppk60(Id2))
X=[xext Ppk60];
Id=kmeans(X,clusterN);
Id1=find(Id==1);
Id2=find(Id==2);
mean(Ppk60(Id2))
std(Ppk60(Id2))
mean(Ppk60(Id1))
std(Ppk60(Id1))
confident
doc confident interval
Z=zscore(Ppk60)
Ppk60
PP=find(~isnan(Ppk60));
PP=Ppk60(PP);
Pz=zscore(PP)
ScatterPlot_widthDistribution
Cellrecord{fn}
ScatterPlot_widthDistribution
close all;Id1=find(Ppk60>=0.02&Ag<=5&Ag>=3);
Id2=find(Ppk60>=0.02&Ag<=9&Ag>=6);
Id3=find(Ppk60>=0.02&Ag<=15&Ag>=10);
P1=[length(Id1)/length(find(Ag<=5&Ag>=3)) length(Id2)/length(find(Ag<=9&Ag>=6)) length(Id3)/length(find(Ag<=15&Ag>=10))]
size(find(Id==1))
size(find(Id==2))
hist(Bnd(:,4))
hist(abs(Bnd(:,4)))
hist(abs(Bnd(Id1,4)))
hist(abs(Bnd(Id2,4)))
hist(abs(Bnd(Id3,4)))
Id1=find(Ag>=3&Ag<=5)
hist(abs(Bnd(Id1,4)))
hold on
Id1=find(Ag>=6&Ag<=8)
hist(abs(Bnd(Id2,4)),'r')
hist(abs(Bnd(Id2,4)))
hold on
Id2=find(Ag>=6&Ag<=8)
hist(abs(Bnd(Id2,4)))
Id3=find(Ag>=9&Ag<=15)
hist(abs(Bnd(Id3,4)))
figure(1)
subplot(3,1,1)
hist(abs(Bnd(Id1,4)))
Id1=find(Ag>=3&Ag<=5)
figure(1)
subplot(3,1,1)
hist(abs(Bnd(Id1,4)))
subplot(2,1,1)
hist(abs(Bnd(Id2,4)))
figure(1)
subplot(3,1,1)
hist(abs(Bnd(Id1,4)))
subplot(3,1,2)
hist(abs(Bnd(Id2,4)))
subplot(3,1,3)
hist(abs(Bnd(Id3,4)))
hist(abs(Bnd(Id3,5)./(Bnd(Id3,5)-Bnd(Id3,4))))
figure(1)
subplot(3,1,1)
hist(abs(Bnd(Id3,5)./(Bnd(Id3,5)-Bnd(Id3,4))))
figure(1)
subplot(3,1,1)
hist(abs(Bnd(Id1,5)./(Bnd(Id1,5)-Bnd(Id1,4))))
subplot(3,1,2)
hist(abs(Bnd(Id2,5)./(Bnd(Id2,5)-Bnd(Id2,4))))
subplot(3,1,3)
hist(abs(Bnd(Id3,5)./(Bnd(Id3,5)-Bnd(Id3,4))))
MB=[mean(abs(Bnd(Id1,5)./(Bnd(Id1,5)-Bnd(Id1,4)))) mean(abs(Bnd(Id2,5)./(Bnd(Id2,5)-Bnd(Id2,4)))) mean(abs(Bnd(Id3,5)./(Bnd(Id3,5)-Bnd(Id3,4))))]
StB=[std(abs(Bnd(Id1,5)./(Bnd(Id1,5)-Bnd(Id1,4)))) std(abs(Bnd(Id2,5)./(Bnd(Id2,5)-Bnd(Id2,4)))) std(abs(Bnd(Id3,5)./(Bnd(Id3,5)-Bnd(Id3,4))))]
errorbar(MB,StB)
errorbar(MB,StB,'k*')
bar([1,2,3],MB)
errorbar(MB,StB,'k*')
hold on
bar([1,2,3],MB)
close all;Id1=find(Ppk60>=0.02&Ag<=6&Ag>=3);
Id2=find(Ppk60>=0.02&Ag<=9&Ag>=7);
Id3=find(Ppk60>=0.02&Ag<=15&Ag>=10);
P1=[length(Id1)/length(find(Ag<=6&Ag>=3)) length(Id2)/length(find(Ag<=9&Ag>=7)) length(Id3)/length(find(Ag<=15&Ag>=10))]
size(find(Id==1))
size(find(Id==2))
figure(100);plot([1,2,3],P1,'r*-','LineWidth',1)
hold on
plot([1,2,3],1-P1,'k*-','LineWidth',1)
close all;Id1=find(Ppk60>=0.05&Ag<=6&Ag>=3);
Id2=find(Ppk60>=0.05&Ag<=9&Ag>=7);
Id3=find(Ppk60>=0.05&Ag<=15&Ag>=10);
P1=[length(Id1)/length(find(Ag<=6&Ag>=3)) length(Id2)/length(find(Ag<=9&Ag>=7)) length(Id3)/length(find(Ag<=15&Ag>=10))]
size(find(Id==1))
size(find(Id==2))
figure(100);plot([1,2,3],P1,'r*-','LineWidth',1)
hold on
plot([1,2,3],1-P1,'k*-','LineWidth',1)
Id1=find(Id==1&Ag<=6&Ag>=3);
Id2=find(Id==1&Ag<=9&Ag>=7);
Id3=find(Id==1&Ag<=15&Ag>=10);
P1=[length(Id1)/length(find(Ag<=6&Ag>=3)) length(Id2)/length(find(Ag<=9&Ag>=7)) length(Id3)/length(find(Ag<=15&Ag>=10))]
size(find(Id==1))
size(find(Id==2))
figure(100);plot([1,2,3],P1,'r*--','LineWidth',1)
hold on
plot([1,2,3],1-P1,'k*-','LineWidth',1)
d1=find(Id==1&Ag<=6&Ag>=3);
Id2=find(Id==1&Ag<=10&Ag>=7);
Id3=find(Id==1&Ag<=15&Ag>=11);
P1=[length(Id1)/length(find(Ag<=6&Ag>=3)) length(Id2)/length(find(Ag<=10&Ag>=7)) length(Id3)/length(find(Ag<=15&Ag>=11))]
size(find(Id==1))
size(find(Id==2))
figure(100);plot([1,2,3],P1,'r*--','LineWidth',1)
hold on
plot([1,2,3],1-P1,'k*-','LineWidth',1)
ScatterPlot_widthDistribution
close all
hist(abs(Bnd(Id1,5)./(Bnd(Id1,5)-Bnd(Id1,4))))
Subplate_relativeposition=abs(Bnd(:,5)./(Bnd(:,5)-Bnd(:,4)));
scatterhist(Subplate_relativeposition, Ppk60,'Group',Id,'Location','SouthEast',...
'Direction','out','Color','kbrm','LineStyle',{'-','-.',':','-'},...
'LineWidth',[2,2,2,2],'Marker','+od*','MarkerSize',[4,5,6,7]);
scatterhist(Subplate_relativeposition, Ppk60)
doc regression
hold on
plotregression(Subplate_relativeposition, Ppk60)
doc regression
regression(Subplate_relativeposition, Ppk60)
regression(Subplate_relativeposition', Ppk60')
confusion(Subplate_relativeposition', Ppk60')
figure(2)
plotconfusion(Subplate_relativeposition, Ppk60)
figure(2)
plotconfusion(Subplate_relativeposition', Ppk60')
Id1=find(Ag<=6&Ag>=3);
plotregression(Subplate_relativeposition(Id1)', Ppk60(Id1)')
plot(Ag,Subplate_relativeposition,'o')
figure(3)
plot(Ag,Subplate_relativeposition,'o')
Id2=find(Ag<=9&Ag>=7);
figure(4);plotregression(Subplate_relativeposition(Id2)', Ppk60(Id2)')
Id3=find(Ag<=15&Ag>=10);
figure(5);plotregression(Subplate_relativeposition(Id3)', Ppk60(Id3)')
plotconfusion(Subplate_relativeposition', Ppksi')
plotregression(Subplate_relativeposition', Ppksi')
X=[xext Ppk60];
Id=kmeans(X,clusterN);
hist(Ppksi(Id==1),0:0.05:1)
figure(5);
hist(Ppksi(Id==1),0:0.05:1)
figure(6)
hist(Ppksi(Id==2),0:0.05:1)
size(find(Id==1))
hist(Ppksi(Id==2),0:0.01:1)
hist(Ppksi(Id==2),0:0.02:1)
hist(Ppksi(Id==2),0:0.001:1)
hist(Ppksi(Id==2),0:0.005:1)
hist(Ppksi(Id==2),0:0.008:1)
hist(Ppksi(Id==2),0:0.02:1)
hist(Ppksi(Id==1),0:0.02:1)
hist(Ppksi(Id==1),0:0.01:1)
figure(6)
figure(7)
hist(Ppksi(Id==2),0:0.01:1)
X=[xext Ppksi];
Id=kmeans(X,clusterN);
scatterhist(xext, Ppksi,'Group',Id,'Location','SouthEast',...
'Direction','out','Color','kbrm','LineStyle',{'-','-.',':','-'},...
'LineWidth',[2,2,2,2],'Marker','+od*','MarkerSize',[4,5,6,7]);
Id=kmeans(X,3);
scatterhist(xext, Ppksi,'Group',Id,'Location','SouthEast',...
'Direction','out','Color','kbrm','LineStyle',{'-','-.',':','-'},...
'LineWidth',[2,2,2,2],'Marker','+od*','MarkerSize',[4,5,6,7]);
ScatterPlot_widthDistribution
close all
ScatterPlot_widthDistribution
nanmean(Ppksi(Id==1))
nanmean(Ppksi(Id==2))
nanmean(Ppksi(Id==1))
Id1=find(Id==1);
nanmean(Ppeaklayersi(Id1,1);)
nanmean(Ppeaklayersi(Id1,1))
nanmean(Ppeaklayersi(find(Id==2),1))
nanmean(Ppeaklayersi(find(Id==3),1))
nanmean(Ppeaklayersi(find(Id==3),3))
nanmean(Ppeaklayersi(find(Id==2),3))
nanmean(Ppeaklayersi(find(Id==1),3))
hist(Bnd(:,3))
figure(100)
hist(Bnd(:,3))
hist(Bnd(:,2))
close all
test
ScatterPlot_widthDistribution
doc corr
%-- 12/3/13, 1:45 PM --%
doc corr
doc image
ScatterPlot_widthDistribution
%-- 12/5/13, 3:30 PM --%
ScatterPlot_widthDistribution
avg60=avg.*256;
avg60(find(avg60>=256))=255;
image(avg60)
image(avg)
C=[avg60;zeros(size(avg60));zeros(size(avg60))]
image(C)
C=[avg60 zeros(size(avg60)) zeros(size(avg60))]
image(C)
doc 3D matrix
doc remap
doc repmat
repmat([1 2; 3 4],[2 3])
repmat([1 2; 3 4],[1 1 3])
C=repmat(avg60,[1 1 3])
C(:,:,2)=0;
image(C)
C=repmat(avg,[1 1 3])
C(:,:,2)=0;
image(C)
avg60=avg;
avg=gridavg2(Xavg,Yavg,Zden50,xtick,ytick,false,Ncell,0);
avg50=avg;
C(:,:3)=avg50
C(:,:,3)=avg50
image(C)
avg=gridavg2(Xavg,Yavg,Zden60,xtick,ytick,false,Ncell,0);
avg60=gridavg2(Xavg,Yavg,Zden60,xtick,ytick,false,Ncell,0);
Ind2=find(avg60<=density); avg60(Ind2)=0;
avg50=gridavg2(Xavg,Yavg,Zden50,xtick,ytick,false,Ncell,0);
Ind2=find(avg50<=density); avg50(Ind2)=0;
C(:,:,1)=avg60;
C(:,:,3)=avg60;
C(:,:,3)=avg50;
image(C)
avg60=gridavg2(Xavg,Yavg,Zden60,xtick,ytick,false,Ncell,0);
Ind2=find(avg60<=density); avg60(Ind2)=0;
avg50=gridavg2(Xavg,Yavg,Zden50,xtick,ytick,false,Ncell,0);
Ind2=find(avg50<=density); avg50(Ind2)=0;
C(:,:,1)=avg60;
C(:,:,2)=0;
C(:,:,3)=avg50;
image(C)
close all
image(C)
C(:,:,1)=avg50;
C(:,:,2)=0;
C(:,:,3)=avg60;
image(C)
C(:,:,1)=avg60;
C(:,:,2)=0;
C(:,:,3)=avg50;
image(C)
ScatterPlot_widthDistribution
avg60=gridavg2(Xavg,Yavg,Zden60,xtick,ytick,false,Ncell,0);
Ind2=find(avg60<=density); avg60(Ind2)=0;
avg50=gridavg2(Xavg,Yavg,Zden50,xtick,ytick,false,Ncell,0);
Ind2=find(avg50<=density); avg50(Ind2)=0;
C(:,:,1)=avg60;
C(:,:,2)=0;
C(:,:,3)=avg50;
image(C)
avg60=gridavg2(Xavg,Yavg,Zden60,xtick,ytick,false,Ncell,0);
Ind2=find(avg60<=density); avg60(Ind2)=0;
avg50=gridavg2(Xavg,Yavg,Zden50,xtick,ytick,false,Ncell,0);
Ind2=find(avg50<=density); avg50(Ind2)=0;
C(:,:,1)=avg60;
C(:,:,2)=0;
C(:,:,3)=avg50;
image(C)
close all
image(C)
ScatterPlot_widthDistribution
avg60=gridavg2(Xavg,Yavg,Zden60,xtick,ytick,false,Ncell,0);
Ind2=find(avg60<=density); avg60(Ind2)=0;
avg50=gridavg2(Xavg,Yavg,Zden50,xtick,ytick,false,Ncell,0);
Ind2=find(avg50<=density); avg50(Ind2)=0;
C(:,:,1)=avg60;
C(:,:,2)=0;
C(:,:,3)=avg50;
image(C)
close all
image(C)
avg60=gridavg2(Xavg,Yavg,Zden60,xtick,ytick,false,Ncell,0);
Ind2=find(avg60<=density); avg60(Ind2)=0;
avg50=gridavg2(Xavg,Yavg,Zden50,xtick,ytick,false,Ncell,0);
Ind2=find(avg50<=density); avg50(Ind2)=0;
C(:,:,1)=avg60;
C(:,:,2)=0;
C(:,:,3)=avg50;
image(C)
close all
image(C)
avg60=gridavg2(Xavg,Yavg,Zden60,xtick,ytick,false,Ncell,0);
Ind2=find(avg60<=density); avg60(Ind2)=0;
avg50=gridavg2(Xavg,Yavg,Zden50,xtick,ytick,false,Ncell,0);
Ind2=find(avg50<=density); avg50(Ind2)=0;
C(:,:,1)=avg60;
C(:,:,2)=0;
C(:,:,3)=avg50;
image(C)
close all
image(C)
test
Meanvalue_Gui_2
Cellrecord{fn}
isfield(Cellrecord{fn},'meanflagPtxInd60')&&isfield(Cellrecord{fn},'meanflagPtxInd50')
a(1).name
a(2).name
guide
%-- 12/5/13, 7:31 PM --%
test
guide
Cellrecord{fn}
jj=find(xx>mxx*Perc);
jj
avg
Cellsum
guide
Cellsum.Bnd
ScatterPlot_widthDistribution
Cellrecord{fn}.DistWidthChargeTp2Hp60
ScatterPlot_widthDistribution
close all
ScatterPlot_widthDistribution
close all
ScatterPlot_widthDistribution
StimcorrX50
ChgCell50
k60
k50
size(ChgCell50{1})
ChgCell50{1}
ChgCell50{2}
ScatterPlot_widthDistribution
ChgCell50{2}
StimcorrX50
k50
k60
MapMatrix=[];
for cn=1:k60
xlim = [-1000-dx-0.1 400];
ylim = [-650-dy-0.1 500];
xtick = xlim(1):dx:xlim(2);
ytick = ylim(1):dy:ylim(2);
xi = xtick(1:end-1)+dx/2;
yi = ytick(1:end-1)+dy/2;
Xavg=StimcorrX60{cn};
Yavg=StimcorrY60{cn};
[cnt_exc_dir ind] =griddensity2(Xavg,Yavg,[],[],xtick,ytick,N,dens);
sel = (cnt_exc_dir<=dens);
isel = find(sel);
selrm = [];
for i=1:length(isel)
sel2 = find(ind==isel(i)); selrm = [selrm sel2];
end
Xavg(selrm) = [];
end
xlim = [-1000-dx-0.1 400];
ylim = [-650-dy-0.1 500];
xtick = xlim(1):dx:xlim(2);
ytick = ylim(1):dy:ylim(2);
xi = xtick(1:end-1)+dx/2;
yi = ytick(1:end-1)+dy/2;
MapMatrix=[];
for cn=1:k60
Xavg=StimcorrX60{cn};
Yavg=StimcorrY60{cn};
Z=ChgCell60{cn};
avg = gridavg2(Xavg,Yavg,Z,xtick,ytick,false,0,0);
Vc=[Vc avg(:)];
end
Vc=[];
a=[1 1; 1 1]
a(:)
size(a(:))
Vc=[];
xlim = [-1000-dx-0.1 400];
ylim = [-650-dy-0.1 500];
xtick = xlim(1):dx:xlim(2);
ytick = ylim(1):dy:ylim(2);
xi = xtick(1:end-1)+dx/2;
yi = ytick(1:end-1)+dy/2;
MapMatrix=[];
for cn=1:k60
Xavg=StimcorrX60{cn};
Yavg=StimcorrY60{cn};
Z=ChgCell60{cn};
avg = gridavg2(Xavg,Yavg,Z,xtick,ytick,false,0,0);
Vc=[Vc avg(:)];
end
size(Vc)
[coeff,score,latent]=pca(Vc);
score
[wcoeff,~,latent,~,explained] = pca(Vc,...
'VariableWeights','variance')
[coeff,score,latent]=pca(Vc);
M=coeff(Vc);
M=coeff(:,1);
size(M)
Vr=latent(:,1)
size(Vr)
size(latent)
size(Vc)
size(score)
[coeff,score,latent]=pca(Vc');
size(coeff)
latent
test_toolbox
help COMPUTE_MAPPPING
help compute_mapping
reduce=compute_mapping(Vc,'PCA')
size(reduce)
size(xi)
size(yi)
M=remat(reduce(:,1),[29, 24])
remap
doc remap
doc remat
M=repmat(reduce(:,1),[29, 24])
image(M)
doc repmat
M=reshape(reduce(:,1),29,24)
image(M)
reduce(:,1)
M=reshape(reduce(:,2),29,24)
imagesc(M)
M=reshape(reduce(:,1),29,24)
imagesc(M)
[coeff,score,latent]=pca(Vc');
size(Vc)
close all
[mappedA, mapping] = compute_mapping(Vc, 'PCA');
recX = reconstruct_data(mappedA, mapping)
imagesc(recX)
imagesc(Vc)
Vc
ScatterPlot_widthDistribution
Vc=[];
xlim = [-1000-dx-0.1 400];
ylim = [-650-dy-0.1 500];
xtick = xlim(1):dx:xlim(2);
ytick = ylim(1):dy:ylim(2);
xi = xtick(1:end-1)+dx/2;
yi = ytick(1:end-1)+dy/2;
MapMatrix=[];
for cn=1:k60
Xavg=StimcorrX60{cn};
Yavg=StimcorrY60{cn};
Z=ChgCell60{cn};
avg = gridavg2(Xavg,Yavg,Z,xtick,ytick,false,0,0);
Vc=[Vc avg(:)];
end
size(Vc)
[mappedA, mapping] = compute_mapping(Vc, 'PCA');
recX = reconstruct_data(mappedA, mapping)
imagesc(recX)
imagesc(Vc)
figure(2)
imagesc(recX)
[mappedA, mapping] = compute_mapping(Vc', 'PCA');
recX = reconstruct_data(mappedA, mapping)
imagesc(recX)
imagesc(Vc)
Vc
imagesc(Vc)
imagesc(recX)
Vc
[coeff,score,latent]=pca(Vc');
size(Vc)
[coeff,score,latent]=pca(Vc');
[coeff,score]=pca(Vc');
size(coeff)
size(score)
score
score.mean
imagesc(score.mean)
imagesc(Vc)
imagesc(Vc')
score
score.lambda
[coeff,score]=pca(Vc');
score
???princomp
help princomp
[coeff,score,latent]=princomp(Vc');
size(coeff)
latent
score
image(coeff(1,:))
image(reshape(coeff(1,:),29,24))
imagesc(reshape(coeff(1,:),29,24))
imagesc(reshape(coeff(1,:)+mean(Vc'),29,24))
imagesc(reshape(coeff(2,:)+mean(Vc'),29,24))
imagesc(reshape(coeff(1,:),29,24))
imagesc(reshape(coeff(1,:)+mean(Vc'),24,29))
imagesc(reshape(coeff(2,:)+mean(Vc'),24,29))
imagesc(reshape(coeff(3,:)+mean(Vc'),24,29))
imagesc(reshape(coeff(4,:)+mean(Vc'),24,29))
imagesc(reshape(coeff(5,:)+mean(Vc'),24,29))
imagesc(reshape(coeff(6,:)+mean(Vc'),24,29))
imagesc(reshape(coeff(7,:)+mean(Vc'),24,29))
xlim = [-1000-dx-0.1 400];
ylim = [-650-dy-0.1 500];
xtick = xlim(1):dx:xlim(2);
ytick = ylim(1):dy:ylim(2);
xi = xtick(1:end-1)+dx/2;
yi = ytick(1:end-1)+dy/2;
imagesc(xi,yi,reshape(coeff(1,:)+mean(Vc'),24,29)); axis image;axis ij;
imagesc(xi,yi,reshape(coeff(2,:)+mean(Vc'),24,29)); axis image;axis ij;
imagesc(xi,yi,reshape(coeff(3,:)+mean(Vc'),24,29)); axis image;axis ij;
imagesc(xi,yi,reshape(coeff(4,:)+mean(Vc'),24,29)); axis image;axis ij;
imagesc(xi,yi,reshape(coeff(6,:)+mean(Vc'),24,29)); axis image;axis ij;
imagesc(xi,yi,reshape(coeff(7,:)+mean(Vc'),24,29)); axis image;axis ij;
imagesc(xi,yi,reshape(coeff(8,:)+mean(Vc'),24,29)); axis image;axis ij;
imagesc(xi,yi,reshape(coeff(9,:)+mean(Vc'),24,29)); axis image;axis ij;
imagesc(xi,yi,reshape(coeff(10,:)+mean(Vc'),24,29)); axis image;axis ij;
imagesc(xi,yi,reshape(coeff(11,:)+mean(Vc'),24,29)); axis image;axis ij;
imagesc(xi,yi,reshape(coeff(12,:)+mean(Vc'),24,29)); axis image;axis ij;
imagesc(xi,yi,reshape(coeff(13,:)+mean(Vc'),24,29)); axis image;axis ij;
imagesc(xi,yi,reshape(coeff(14,:)+mean(Vc'),24,29)); axis image;axis ij;
imagesc(xi,yi,reshape(coeff(15,:)+mean(Vc'),24,29)); axis image;axis ij;
size(score)
imagesc(xi,yi,reshape(score(1,:)+mean(Vc'),24,29)); axis image;axis ij;
imagesc(xi,yi,reshape(score(2,:)+mean(Vc'),24,29)); axis image;axis ij;
imagesc(xi,yi,reshape(score(3,:)+mean(Vc'),24,29)); axis image;axis ij;
imagesc(xi,yi,reshape(score(4,:)+mean(Vc'),24,29)); axis image;axis ij;
imagesc(xi,yi,reshape(mean(Vc'),24,29)); axis image;axis ij;
imagesc(xi,yi,reshape(coeff(15,:),24,29)); axis image;axis ij;
imagesc(xi,yi,reshape(coeff(1,:),24,29)); axis image;axis ij;
imagesc(xi,yi,reshape(coeff(2,:),24,29)); axis image;axis ij;
imagesc(xi,yi,reshape(coeff(3,:),24,29)); axis image;axis ij;
imagesc(xi,yi,reshape(coeff(4,:),24,29)); axis image;axis ij;
imagesc(xi,yi,reshape(coeff(5,:),24,29)); axis image;axis ij;
mu=mean(Vc');
[w pc ev]=princomp(Vc');
xhat=bsxfun(@minus,Vc',mu);
XX=mu+pc*w';
size(w)
size(pc)
size(mu)
XX=mu+pc*w'(1,:);
XX=mu+pc(1,:)*w'(1,:);
XX=mu+pc(1,:)*w';
imagesc(xi,yi,reshape(XX,24,29)); axis image;axis ij;
XX=mu+pc(2,:)*w';
imagesc(xi,yi,reshape(XX,24,29)); axis image;axis ij;
XX=mu+pc(3,:)*w';
imagesc(xi,yi,reshape(XX,24,29)); axis image;axis ij;
XX=mu+pc(4,:)*w';
imagesc(xi,yi,reshape(XX,24,29)); axis image;axis ij;
XX=mu+pc(5,:)*w';
imagesc(xi,yi,reshape(XX,24,29)); axis image;axis ij;
XX=mu+pc(6,:)*w';
imagesc(xi,yi,reshape(XX,24,29)); axis image;axis ij;
XX=mu+pc(7,:)*w';
imagesc(xi,yi,reshape(XX,24,29)); axis image;axis ij;
XX=mu+pc(8,:)*w';
imagesc(xi,yi,reshape(XX,24,29)); axis image;axis ij;
XX=mu+pc(9,:)*w';
latent
XX=mu+pc(1,:)*w';
imagesc(xi,yi,reshape(XX,24,29)); axis image;axis ij;
XX=mu+pc(2,:)*w';
imagesc(xi,yi,reshape(XX,24,29)); axis image;axis ij;
xapprox=pc(:,1:2)*w(:,1:2)';
size(xapprox)
figure(1)
subplot(2,1,1)
imagesc(Vc')
xapprox=bsxfun(@plus,mu,xapprox)
subplot(2,1,2)
imagesc(xapprox)
xapprox=pc(:,1:4)*w(:,1:4)';
xapprox=bsxfun(@plus,mu,xapprox)
subplot(2,1,1)
imagesc(Vc')
subplot(2,1,2)
imagesc(xapprox)
XX=mu+pc(1,:)*w';
imagesc(xi,yi,reshape(XX,24,29)); axis image;axis ij;
figure(2)
XX=mu+pc(2,:)*w';
imagesc(xi,yi,reshape(w(:,1)'+mu,24,29)); axis image;axis ij;
imagesc(xi,yi,reshape(w(:,2)'+mu,24,29)); axis image;axis ij;
imagesc(xi,yi,reshape(w(:,3)'+mu,24,29)); axis image;axis ij;
imagesc(xi,yi,reshape(w(:,3)'24,29)); axis image;axis ij;
imagesc(xi,yi,reshape(w(:,1)',24,29)); axis image;axis ij;
imagesc(xi,yi,reshape(w(:,2)',24,29)); axis image;axis ij;
imagesc(xi,yi,reshape(w(:,3)',24,29)); axis image;axis ij;
imagesc(xi,yi,reshape(w(:,4)',24,29)); axis image;axis ij;
imagesc(xi,yi,reshape(pc(:,1)',24,29)); axis image;axis ij;
imagesc(xi,yi,reshape(pc(1,:),24,29)); axis image;axis ij;
imagesc(xi,yi,reshape(pc(2,:),24,29)); axis image;axis ij;
imagesc(xi,yi,reshape(pc(3,:),24,29)); axis image;axis ij;
imagesc(xi,yi,reshape(pc(4,:),24,29)); axis image;axis ij;
imagesc(xi,yi,reshape(w(:,5)',24,29)); axis image;axis ij;
imagesc(xi,yi,reshape(w(:,6)',24,29)); axis image;axis ij;
imagesc(xi,yi,reshape(w(:,7)',24,29)); axis image;axis ij;
size(xapprox)
imagesc(xi,yi,reshape(xapprox(:,1)',24,29)); axis image;axis ij;
imagesc(xi,yi,reshape(xapprox(1,:),24,29)); axis image;axis ij;
imagesc(xi,yi,reshape(Vc(:,1)',24,29)); axis image;axis ij;
imagesc(xi,yi,reshape(xapprox(:,1)',24,29)); axis image;axis ij;
figure(1)
subplote(2,1,1)
subplot(2,1,1)
imagesc(xi,yi,reshape(Vc(:,1)',24,29)); axis image;axis ij;
subplot(2,1,2)
imagesc(xi,yi,reshape(xapprox(1,:),24,29)); axis image;axis ij;
for i=1:105
figure(i)
> subplot(2,1,1)
imagesc(xi,yi,reshape(Vc(:,1)',24,29)); axis image;axis ij;
subplot(2,1,2)
imagesc(xi,yi,reshape(xapprox(1,:),24,29)); axis image;axis ij;
for i=1:105
figure(i)
subplot(2,1,1)
imagesc(xi,yi,reshape(Vc(:,i)',24,29)); axis image;axis ij;
subplot(2,1,2)
imagesc(xi,yi,reshape(xapprox(i,:),24,29)); axis image;axis ij;
end
xapprox=pc(:,1:6)*w(:,1:6)';
xapprox=bsxfun(@plus,mu,xapprox)
for i=1:105
figure(i)
subplot(2,1,1)
imagesc(xi,yi,reshape(Vc(:,i)',24,29)); axis image;axis ij;
subplot(2,1,2)
imagesc(xi,yi,reshape(xapprox(i,:),24,29)); axis image;axis ij;
end
close all
xapprox=pc(:,:)*w(:,:)';
xapprox=bsxfun(@plus,mu,xapprox)
for i=1:105
figure(i)
subplot(2,1,1)
imagesc(xi,yi,reshape(Vc(:,i)',24,29)); axis image;axis ij;
subplot(2,1,2)
imagesc(xi,yi,reshape(xapprox(i,:),24,29)); axis image;axis ij;
end
xapprox=pc(:,1:10)*w(:,1:10)';
xapprox=bsxfun(@plus,mu,xapprox)
for i=1:105
figure(i)
subplot(2,1,1)
imagesc(xi,yi,reshape(Vc(:,i)',24,29)); axis image;axis ij;
subplot(2,1,2)
imagesc(xi,yi,reshape(xapprox(i,:),24,29)); axis image;axis ij;
end
close all
doc imagesc
ScatterPlot_widthDistribution
C
find(C>1)
find(C<0)
C(find(C<0))
ScatterPlot_widthDistribution
close all
ScatterPlot_widthDistribution
%-- 12/18/13, 1:33 PM --%
ScatterPlot_widthDistribution
C
find(C>1)
ScatterPlot_widthDistribution
Ag
ScatterPlot_widthDistribution
test
ScatterPlot_widthDistribution
close all
ScatterPlot_widthDistribution
doc kmean
ScatterPlot_widthDistribution
xx=[0.5217 0.3913 0.0435 0.0435;0.3542 0.2083 0.2917 0.1458;0.3684 0.1316 0.3158 0.1842]
hist(xx)
bar(xx)
ScatterPlot_widthDistribution
xx=[0.6087 0.3478 0.0435; 0.4375 0.3542 0.2083;0.3684 0.3158 0.3158]
bar(xx)
yy(:,1)=xx;
yy(:,1)=xx(:,1);
yy(:,2)=xx(:,3);
yy(:,3)=xx(:,2);
bar(yy)
close all
bar(yy)
ScatterPlot_widthDistribution
close all
ScatterPlot_widthDistribution
xx=[0.6522 0.3478;0.6255 0.3750;0.6842 0.3158]
bar(xx)
plot(xx,'*-')
xx=[0.6087 0.3478 0.0435; 0.4375 0.3542 0.2083;0.3684 0.3158 0.3158]
plot(xx,'*-')
plot(yy,'*-')
close all
plot(yy,'*-')
xx=[0.6522 0.3478;0.6255 0.3750;0.6842 0.3158]
plot(xx,'*-')
ScatterPlot_widthDistribution
xx=[0.5217 0.3913 0.0435 0.0435;0.3542 0.2083 0.2917 0.1458;0.3684 0.1316 0.3158 0.1842]
plot(xx,'*-')
bar(xx)
ScatterPlot_widthDistribution
close all
ScatterPlot_widthDistribution
close all
ScatterPlot_widthDistribution
xlim([-475 475])
doc xlim
xlim([-475 475])
23/109
15/109
28/109
43/109
43+28+15+23
28/109
49/109
ScatterPlot_widthDistribution
close all
ScatterPlot_widthDistribution
close all
ScatterPlot_widthDistribution
close all
ScatterPlot_widthDistribution
xx=[ 0.5217    0.3913    0.0435    0.0435;   0.3542    0.2917   0.2083  0.1458]
ScatterPlot_widthDistribution
xx=[xx;0.3684   0.1316  0.3158    0.1842]
plot(xx,'*-')
ScatterPlot_widthDistribution
close all
ScatterPlot_widthDistribution
close all
ScatterPlot_widthDistribution
%-- 12/23/13, 12:33 PM --%
guide
close all
test
guifr
guide
62*log(15/137)
62*log(137/15)
doc log
62*log10(137/15)
figure(100)
plot(resampdata(:,1))
max(max(resampdata))-min(min(resampdata))
max(resampdata)
min(resampdata)
max(resampdata)
max(max(resampdata))-min(min(resampdata))
figure(100)
plot(resampdata(:,1))
figure(100)
plot(resampdata(:,1))
plot(resampdata(:,2))
figure(100)
p
plot(resampdata(:,2))
plot(resampdata - repmat(min(resampdata),[rsrate-9 1]))
size(resampdata)
figure(100)
plot(resampdata - repmat(min(resampdata),[rsrate-9 1]))
plot((resampdata - repmat(min(resampdata),[rsrate-9 1]))./(max(max(resampdata))-min(min(resampdata))))
figure(100)
plot(resampdata(:,1))
plot(resampdata(:,:))
alignedhistWithMaps(/Volumes/disk1/silent_122313/eventwindow50')
alignedhistWithMaps('/Volumes/disk1/silent_122313/eventwindow50')
test
figure(100)
plot(resampdata)
figure(100)
plot(resampdata)
ScatterPlot_widthDistribution
%-- 12/31/13, 10:45 AM --%
ScatterPlot_widthDistribution
close all
%-- 12/31/13, 11:17 AM --%
ScatterPlot_widthDistribution
%-- 12/31/13, 11:23 AM --%
ScatterPlot_widthDistribution
[cid1,cid2]=find(avg2_60>0|avg2_50>0);
cid1
size(cid1)
size(cid2)
%-- 12/31/13, 11:32 AM --%
ScatterPlot_widthDistribution
size(cid2)
size(cid1)
size(C)
doc find
ScatterPlot_widthDistribution
size(avg2_60)
ScatterPlot_widthDistribution
size(cid1)
size(avg2_60+avg2_50)
doc find
C(cid1,cid2,2)=0;
%-- 12/31/13, 11:46 AM --%
ScatterPlot_widthDistribution
silentsynapses_Gui
ScatterPlot_widthDistribution
C2
ScatterPlot_widthDistribution
size(C2)
size(avg60)
ScatterPlot_widthDistribution
%-- 1/9/14, 10:44 AM --%
load '/Volumes/disk1/silent_122313/eventwindow50/Map_121001_1634.mat'
load '/Volumes/disk1/silent_122313/eventwindow50/LSPS_Map_121001_1634.mat'
cells
events
cells.events
La=[];DU=[];
AREA=[];PA=[];TAU=[];eventcnt=0;Flg=[];
for k=1:cells.header.nPts
%for dd=1:1:length(Indx)
%CL=CC{dd};
header1=cells.header;
stim_start=find(header1.stimulatorOutput(4,:));
%xxx=INMatrix(k,dd);
fevent=find(cells.events.flag{k}(:)>0);
if length(fevent)
if cells.events.flag{k}(fevent(1))==1
ff=NaN;
elseif cells.events.flag{k}(fevent(1))==2
ff=fevent(1);
end
else ff=NaN;
end
if ~isnan(ff)
latency=cells.events.startSamp{k}(ff).*1000./cells.header.sampleRate-stim_start(1).*1000./header1.sampleRate;
dur=cells.events.duration{k}(ff).*1000./CL.cells.header.sampleRate;
area=(cells.events.area{k}(ff));
peak=(cells.events.peakAmp{k}(ff));
tau=(cells.events.tau{k}(ff));
flag_latency=1;
%            eventcnt=eventcnt+1;
else
latency=NaN;
dur=NaN;
area=NaN;
peak=NaN;
tau=NaN;
flag_latency=0;
end
La=[La,latency];
DU=[DU dur];
AREA=[AREA area];
PA=[PA peak];
TAU=[TAU tau];
Flg=[Flg flag_latency];
end
La=[];DU=[];
AREA=[];PA=[];TAU=[];eventcnt=0;Flg=[];
for k=1:cells.header.nPts
%for dd=1:1:length(Indx)
%CL=CC{dd};
header1=cells.header;
stim_start=find(header1.stimulatorOutput(4,:));
%xxx=INMatrix(k,dd);
fevent=find(cells.events.flag{k}(:)>0);
if length(fevent)
if cells.events.flag{k}(fevent(1))==1
ff=NaN;
elseif cells.events.flag{k}(fevent(1))==2
ff=fevent(1);
end
else ff=NaN;
end
if ~isnan(ff)
latency=cells.events.startSamp{k}(ff).*1000./cells.header.sampleRate-stim_start(1).*1000./header1.sampleRate;
dur=cells.events.duration{k}(ff).*1000./cells.header.sampleRate;
area=(cells.events.area{k}(ff));
peak=(cells.events.peakAmp{k}(ff));
tau=(cells.events.tau{k}(ff));
flag_latency=1;
%            eventcnt=eventcnt+1;
else
latency=NaN;
dur=NaN;
area=NaN;
peak=NaN;
tau=NaN;
flag_latency=0;
end
La=[La,latency];
DU=[DU dur];
AREA=[AREA area];
PA=[PA peak];
TAU=[TAU tau];
Flg=[Flg flag_latency];
end
nanmean(PA)
Chg=AREA;
PK=PA;
load '/Volumes/disk1/silent_122313/eventwindow50/LSPS_Map_121001_1706.mat'
La=[];DU=[];
AREA=[];PA=[];TAU=[];eventcnt=0;Flg=[];
for k=1:cells.header.nPts
%for dd=1:1:length(Indx)
%CL=CC{dd};
header1=cells.header;
stim_start=find(header1.stimulatorOutput(4,:));
%xxx=INMatrix(k,dd);
fevent=find(cells.events.flag{k}(:)>0);
if length(fevent)
if cells.events.flag{k}(fevent(1))==1
ff=NaN;
elseif cells.events.flag{k}(fevent(1))==2
ff=fevent(1);
end
else ff=NaN;
end
if ~isnan(ff)
latency=cells.events.startSamp{k}(ff).*1000./cells.header.sampleRate-stim_start(1).*1000./header1.sampleRate;
dur=cells.events.duration{k}(ff).*1000./cells.header.sampleRate;
area=(cells.events.area{k}(ff));
peak=(cells.events.peakAmp{k}(ff));
tau=(cells.events.tau{k}(ff));
flag_latency=1;
%            eventcnt=eventcnt+1;
else
latency=NaN;
dur=NaN;
area=NaN;
peak=NaN;
tau=NaN;
flag_latency=0;
end
La=[La,latency];
DU=[DU dur];
AREA=[AREA area];
PA=[PA peak];
TAU=[TAU tau];
Flg=[Flg flag_latency];
end
PKAPV=PA;
ChgAPV=AREA;
Pk60=PK;
Chg60=Chg;
PkAPV60=PKAPV;
ChgAPV60=ChgAPV;
load '/Volumes/disk1/silent_122313/eventwindow50/LSPS_Map_121001_1639.mat'
La=[];DU=[];
AREA=[];PA=[];TAU=[];eventcnt=0;Flg=[];
for k=1:cells.header.nPts
%for dd=1:1:length(Indx)
%CL=CC{dd};
header1=cells.header;
stim_start=find(header1.stimulatorOutput(4,:));
%xxx=INMatrix(k,dd);
fevent=find(cells.events.flag{k}(:)>0);
if length(fevent)
if cells.events.flag{k}(fevent(1))==1
ff=NaN;
elseif cells.events.flag{k}(fevent(1))==2
ff=fevent(1);
end
else ff=NaN;
end
if ~isnan(ff)
latency=cells.events.startSamp{k}(ff).*1000./cells.header.sampleRate-stim_start(1).*1000./header1.sampleRate;
dur=cells.events.duration{k}(ff).*1000./cells.header.sampleRate;
area=(cells.events.area{k}(ff));
peak=(cells.events.peakAmp{k}(ff));
tau=(cells.events.tau{k}(ff));
flag_latency=1;
%            eventcnt=eventcnt+1;
else
latency=NaN;
dur=NaN;
area=NaN;
peak=NaN;
tau=NaN;
flag_latency=0;
end
La=[La,latency];
DU=[DU dur];
AREA=[AREA area];
PA=[PA peak];
TAU=[TAU tau];
Flg=[Flg flag_latency];
end
Pk40=PA;
Chg40=AREA;
load '/Volumes/disk1/silent_122313/eventwindow50/LSPS_Map_121001_1711.mat'
La=[];DU=[];
AREA=[];PA=[];TAU=[];eventcnt=0;Flg=[];
for k=1:cells.header.nPts
%for dd=1:1:length(Indx)
%CL=CC{dd};
header1=cells.header;
stim_start=find(header1.stimulatorOutput(4,:));
%xxx=INMatrix(k,dd);
fevent=find(cells.events.flag{k}(:)>0);
if length(fevent)
if cells.events.flag{k}(fevent(1))==1
ff=NaN;
elseif cells.events.flag{k}(fevent(1))==2
ff=fevent(1);
end
else ff=NaN;
end
if ~isnan(ff)
latency=cells.events.startSamp{k}(ff).*1000./cells.header.sampleRate-stim_start(1).*1000./header1.sampleRate;
dur=cells.events.duration{k}(ff).*1000./cells.header.sampleRate;
area=(cells.events.area{k}(ff));
peak=(cells.events.peakAmp{k}(ff));
tau=(cells.events.tau{k}(ff));
flag_latency=1;
%            eventcnt=eventcnt+1;
else
latency=NaN;
dur=NaN;
area=NaN;
peak=NaN;
tau=NaN;
flag_latency=0;
end
La=[La,latency];
DU=[DU dur];
AREA=[AREA area];
PA=[PA peak];
TAU=[TAU tau];
Flg=[Flg flag_latency];
end
Pk40APV=PA;
Chg40APV=AREA;
mPK60=nanmean(Pk60);
mChg60=nanmean(Chg60);
mPK40=nanmean(Pk40);
mChg40=nanmean(Chg40);
mPK60APV=nanmean(Pk60APV);
mChg60APV=nanmean(Chg60APV);
mPK40APV=nanmean(Pk40APV);
mChg40APV=nanmean(Chg40APV);
mPK60APV=nanmean(PK60APV);
mPK60APV=nanmean(PkAPV60);
mChg60APV=nanmean(ChgAPV60);
mPK40APV=nanmean(Pk40APV);
mChg40APV=nanmean(Chg40APV);
bar([mPK60 mPK60APV mPK40 mPK40APV])
doc bar
bar([mPK60 mPK60APV ; mPK40 mPK40APV])
vPK60=nanstd(Pk60);
vChg60=nanstd(Chg60);
vPK40=nanstd(Pk40);
vChg40=nanstd(Chg40);
vPK60APV=nanstd(PkAPV60);
vChg60APV=nanstd(ChgAPV60);
vPK40APV=nanstd(Pk40APV);
vChg40APV=nanstd(Chg40APV);
error([vPK60 vPK60APV vPK40 VPK40APV])
error([vPK60 vPK60APV vPK40 vPK40APV])
error([mPK60 mPK60APV mPK40 mPK40APV],[vPK60 vPK60APV vPK40 vPK40APV])
doc error
doc errorbar
errorbar([mPK60 mPK60APV mPK40 mPK40APV],[vPK60 vPK60APV vPK40 vPK40APV])
x=[0.8 1.2 1.8 2.2];
errorbar(x,[mPK60 mPK60APV mPK40 mPK40APV],[vPK60 vPK60APV vPK40 vPK40APV])
hold on
bar([mPK60 mPK60APV ; mPK40 mPK40APV])
x=[0.85 1.15 1.85 2.15];
errorbar(x,[mPK60 mPK60APV mPK40 mPK40APV],[vPK60 vPK60APV vPK40 vPK40APV])
hold on
bar([mPK60 mPK60APV ; mPK40 mPK40APV])
errorbar(x,[mPK60 mPK60APV mPK40 mPK40APV],[vPK60 vPK60APV vPK40 vPK40APV],'b*')
hold on
bar([mPK60 mPK60APV ; mPK40 mPK40APV])
x=[0.87 1.13 1.87 2.13];
errorbar(x,[mPK60 mPK60APV mPK40 mPK40APV],[vPK60 vPK60APV vPK40 vPK40APV],'b*')
hold on
bar([mPK60 mPK60APV ; mPK40 mPK40APV])
x=[0.86 1.14 1.86 2.14];
errorbar(x,[mPK60 mPK60APV mPK40 mPK40APV],[vPK60 vPK60APV vPK40 vPK40APV],'b*')
hold on
bar([mPK60 mPK60APV ; mPK40 mPK40APV])
%-- 1/10/14, 1:42 PM --%
load '/Volumes/disk1/silent_122313/eventwindow50/LSPS_Map_121001_1634.mat'
1-NaN
AVGPEAK_121001
hold on
plot(0,0,'o')
hold on
plot(0,0,'o')
RM_vwhCombine_Eular_constant_noise
AVGPEAK_121001
close all
AVGPEAK_121001
hold on
plot(0,0,'wo')
hold on
plot(0,0,'wo')
AVGPEAK_121001
close all
AVGPEAK_121001
close all
RM_vwhCombine_Eular_constant_noise
%-- 1/12/14, 11:02 AM --%
RM_vwhCombine_Eular_constant_noise
log2
log10(2)
hold on
figure(2)
hold on
v=-110:1:20;
v0=-110:1:20;
m0=(1+exp(-(v0+38)./7)).^(-1);
h0=(1+exp((v0+65+6)./6))^(-1);
h0=(1+exp((v0+65+6)./6)).^(-1);
w0=(1+exp(-(v0+48)./6)).^(-1/4);
w04=w0.^4;
m03=m0.^3;
plot(v0,m0,v0,m03,'--',v0,h0,v0,w0,v0,w04,'--')
figure(3)
plot(v0,m0,v0,m03,'--',v0,h0,v0,w0,v0,w04,'--')
test
close(figure(4))
close all
guide
test
close all
ScatterPlot_widthDistribution_tmc
Bnd
Id1
mC
ScatterPlot_widthDistribution_tmc
mean(Bnd(Id1,:))
Id1
Bnd
ScatterPlot_widthDistribution_tmc
ScatterPlot_widthDistribution
ScatterPlot_widthDistribution_tmc
silentsynapses_Gui
ScatterPlot_widthDistribution
firingrate_poisson_vwh_combine
AVGPEAK_121001
hist(Pk40)
hist(Pk40APV)
figure(17)
hist(Pk40APV)
anova1(Pk40, Pk40APV)
anova1(Pk60, Pk60APV)
anova1(Pk60, PkAPV60)
doc ttest
[h,p]=ttest(Pk60-PkAPV60)
[h,p]=ttest(Pk40-PkAPV40)
[h,p]=ttest(Pk40-Pk40APV)
hist(Pk40APV)
lillietest(Pk40APV)
lillietest(Pk40)
signtest(Pk40,Pk40APV)
signrank(Pk40,Pk40APV)
signrank(Pk60,PkAPV60)
lillietest(Pk40-Pk40APV)
ttest(Pk40-Pk40APV,0)
[h,p]=ttest(Pk40-Pk40APV)
test
ScatterPlot_widthDistribution
CorrPtx60
ScatterPlot_widthDistribution
Cellrecord{1}
ScatterPlot_widthDistribution
dir(fullfile(folderMeanpath,'*.mat'))
ScatterPlot_widthDistribution
Cellrecord{1}
CorrPtx60
hist(abs(CorrPtx60))
figure(20)
hist(abs(CorrPtx60))
hist(abs(CorrPtx10))
hist(abs(CorrPtx50))
CorrPtx50
hist(abs(CorrPtx60))
mean(CorrPtx60)
nanmean(CorrPtx60)
nanmean(abs(CorrPtx60))
nanstd(abs(CorrPtx60))
%-- 1/21/14, 10:07 AM --%
firingrate_poisson_vwh_combine
load Poisson_firingrate_vwhCombine_Gmax10_transientinh0_hshif6modifiedNa.mat
figure(13)
plot(Vmean,Fr)
plot(Vmean,Fr.*1000)
figure(1)
plot(RR.*1000,Fr.*1000)
hold on
figure(2)
plot(Vmean,Fr.*1000)
hold on
load Poisson_firingrate_vwhCombine_Gmax10_transientinh10_hshif6modifiedNa.mat
figure(1)
plot(RR.*1000,Fr.*1000,r)
hold on
figure(2)
plot(Vmean,Fr.*1000,r)
hold on
figure(1)
plot(RR.*1000,Fr.*1000,'r')
hold on
figure(2)
plot(Vmean,Fr.*1000,'r')
hold on
load Poisson_firingrate_vwhCombine_Gmax10_transientinh20_hshif6modifiedNa.mat
figure(1)
plot(RR.*1000,Fr.*1000,'g')
hold on
figure(2)
plot(Vmean,Fr.*1000,'g')
hold on
figure(3)
plot(RR,Vmean,'g')
hold on
load Poisson_firingrate_vwhCombine_Gmax10_transientinh10_hshif6modifiedNa.mat
plot(RR,Vmean,'g')
plot(RR,Vmean,'r')
load Poisson_firingrate_vwhCombine_Gmax10_transientinh0_hshif6modifiedNa.mat
plot(RR,Vmean,'b')
test
%-- 1/25/14, 10:28 AM --%
test
ScatterPlot_widthDistribution_tmc
ScatterPlot_widthDistribution_darkexp
mean(Bnd(Id1,:),1)
ScatterPlot_widthDistribution_darkexp
%-- 1/26/14, 6:53 PM --%
test
%-- 1/26/14, 6:58 PM --%
test
LSPSgetevents3
%-- 1/27/14, 12:18 PM --%
test
[]
cells
gui
guide
doc anova
h=norm(100,5)
doc norm
h=rand(100,5)
[p,tb1,stats]=anova1(h)
[c,m] = multcompare(stats)
stats
ScatterPlot_widthDistribution_darkexp
close all
doc kruskalwallis
doc multicomp
x=ones(100,1);
x=[x,x.*2,x.*3];
x=x+rand(100,3)
[p,tbl,stats]=kruskalwallis(x);
multcomp(stats,'estimate','kruskalwallis')
multcompare(stats,'estimate','kruskalwallis')
ScatterPlot_widthDistribution_darkexp
doc line
line([-50 50],[0 0])
figure(100)
line([-50 50],[0 0])
ScatterPlot_widthDistribution_darkexp
doc
doc line
ScatterPlot_widthDistribution_darkexp
close all
test
fid = fopen('/Users/xymeng/Documents/MATLAB/drtoolbox 2/Readme.txt');
fid
fid = importdata('/Users/xymeng/Documents/MATLAB/drtoolbox 2/Readme.txt');
fid
fid{1}
fid.data
fid.textdata
fid = importdata('/Users/xymeng/Documents/MATLAB/drtoolbox 2/Readme.txt',delimiterIn);
fid = importdata('/Users/xymeng/Documents/MATLAB/drtoolbox 2/Readme.txt',1);
fid = importdata('/Users/xymeng/Documents/MATLAB/drtoolbox 2/Readme.txt', ,1);
delimiterIn=' ';
fid = importdata('/Users/xymeng/Documents/MATLAB/drtoolbox 2/Readme.txt',delimiterIn,1);
fid
fid = importdata('/Users/xymeng/Documents/MATLAB/drtoolbox 2/Readme.txt',delimiterIn,50);
fid
fid(1)
fid(2)
for i=1:50
c=fid(i);
strcmp(c,'Version')
end
for i=1:10
c=fid(i);
strcmp(c,'Version')
end
doc strcmp
for i=1:50
c=fid(i);
findstr(c,'Version')
end
doc findstr
string(fid(i))
str(fid(i))
fid(i)
findstr(fid(i),'Matlab')
char(fid(i))
findstr(char(fid(i)),'Matlab')
for i=1:50
c=fid(i);
if findstr(char(c),'Version:')
end
end
for i=1:50
c=fid(i);
if findstr(char(c),'Version:') n=indstr(char(c),'Version:'); t=char(c); v=t(n+2);
end
end
for i=1:50
c=fid(i);
if findstr(char(c),'Version:') n=indstr(char(c),'Version:'); t=char(c); v=t(n+2,n+3);
end
end
for i=1:50
c=fid(i);
if findstr(char(c),'Version:') n=indstr(char(c),'Version:'); v=c(n+2);
end
end
doc indstr
for i=1:50
c=fid(i);
if findstr(char(c),'Version:') n=indstr(char(c),'Version:'); v=t(n+2:n+3);
end
end
char(c)
c(1)
for i=1:50
c=fid(i);
if findstr(char(c),'Version:') n=findstr(char(c),'Version:'); g=regexp(c,'','split'); v=g(n+1)
end
end
g=regexp(c,'','split')
g{1}
g=strread(c,'%s')
doc strread
g=strread(char(c),'%s')
for i=1:50
c=fid(i);
if findstr(char(c),'Version:') n=findstr(char(c),'Version:'); t=strread(char(c),'%s'); v=t(n+2);
end
end
for i=1:50
c=fid(i);
if findstr(char(c),'Version:') n=findstr(char(c),'Version:'); t=strread(char(c),'%s'); v=t(n+1);
end
end
v
fid = importdata('/Users/xymeng/Documents/MATLAB/drtoolbox 2/Readme.txt',delimiterIn,Inf);
doc Fread
fid = fopen('/Users/xymeng/Documents/MATLAB/drtoolbox 2/Readme.txt');
m=fread(fid)
m(1)
fid = fopen('/Users/xymeng/Documents/MATLAB/drtoolbox 2/Readme.txt');
fread(fid)
fid = importdata('/Users/xymeng/Documents/MATLAB/drtoolbox 2/Readme.txt',delimiterIn,100);
n=100;
fid = importdata('/Users/xymeng/Documents/MATLAB/drtoolbox 2/Readme.txt',delimiterIn,n);
for i=1:n
c=fid(i);
nstr=findstr(char(c),'Version:');
if nstr
cstr=strread(char(c),'%s');
vstr=cstr(nstr+1);
end
end
vstr
ReadoutFromTxt_forDan
A=[1 3; 2 4]
dataset2cell(A)
>>
doc char
pd = makedist('Binomial')
pd
test
guide
ScatterPlot_widthDistribution_darkexp
%-- 2/17/14, 10:19 AM --%
Cellresponds_gui
header=cells.header;
[imd imdIdx] = sort(header.ImageData(:));
imd(imdIdx) = floor(linspace(0,256-eps,numel(header.ImageData)));
imd = uint8(reshape(imd, size(header.ImageData)));
imAxes=subplot(1,1,1);
ImageScale=header.ImageScale;
image([-1 1].*ImageScale(1)/2, [-1 1].*ImageScale(2)/2, repmat(imd,[1 1 3]),'Parent',imAxes);
axis xy;
set(gca,'XTick',[-1 -.5 0 .5 1].*ImageScale(1)/2, 'YTick',[-1 -.5 0 .5 1].*ImageScale(2)/2, ...
'DataAspectRatio',[ 1 1 1]);
hold on;
header=cells.header;
[imd imdIdx] = sort(header.ImageData(:));
imd(imdIdx) = floor(linspace(0,256-eps,numel(header.ImageData)));
imd = uint8(reshape(imd, size(header.ImageData)));
imAxes=subplot(1,1,1);
ImageScale=header.ImageScale;
image([-1 1].*ImageScale(2)/2, [-1 1].*ImageScale(1)/2, repmat(imd',[1 1 3]),'Parent',imAxes);
axis xy;
set(gca,'XTick',[-1 -.5 0 .5 1].*ImageScale(2)/2, 'YTick',[-1 -.5 0 .5 1].*ImageScale(1)/2, ...
'DataAspectRatio',[ 1 1 1]);
hold on;
[imd imdIdx] = sort(header.ImageData(:));
imd(imdIdx) = floor(linspace(0,256-eps,numel(header.ImageData)));
imd = uint8(reshape(imd, size(header.ImageData)));
imAxes=subplot(1,1,1);
ImageScale=header.ImageScale;
image([-1 1].*ImageScale(2)/2, [-1 1].*ImageScale(1)/2, repmat(imd',[1 1 3]),'Parent',imAxes);
axis ij;
set(gca,'XTick',[-1 -.5 0 .5 1].*ImageScale(2)/2, 'YTick',[-1 -.5 0 .5 1].*ImageScale(1)/2, ...
'DataAspectRatio',[ 1 1 1]);
hold on;
plot(header.StimCoordinates(2,:),header.StimCoordinates(1,:),'b.')
test
ScatterPlot_widthDistribution_darkexp
close all
ScatterPlot_widthDistribution_darkexp
close all
ScatterPlot_widthDistribution_darkexp
close all
avghistLayers
ScatterPlot_widthDistribution_darkexp
D60x=D60pk(:,9);
D60x=D60Pk(:,9);
hsit(D60x)
hist(D60x)
close all
hist(D60x)
D60x=D60Pk(:,10);
hist(D60x)
D60x=D60Pk(:,11);
D60x=D60Pk(:,9);
D60x_Norm=D60x;
D10x=D10Pk(:,9);
D10x_Norm=D10x;
ScatterPlot_widthDistribution_darkexp
D60x=D60Pk(:,9);
D10x=D10Pk(:,9);
D60x_DE=D60x;
D10x_DE=D10x;
hist(D60x_DE)
kstest(D60x_Norm,D60x_DE)
doc kstest
size(D60x_Norm)
size(D60x_DE)
hist(D60x_DE)
close all
hist(D60x_DE)
figure(2)
hist(D60x_Norm)
kstest(D60x_Norm-D60x_DE)
[h,p]=kstest(D60x_Norm-D60x_DE)
doc boxplot
mean(D60x_Norm)
nanmean(D60x_Norm)
D60x_Norm(find(isnan(D60x_Norm)))=0;
D60x_DE(find(isnan(D60x_DE)))=0;
mean(D60x_Norm)
mean(D60x_DE)
mDE=mean(D60x_DE);
vDE=std(D60x_DE);
D10x_DE(find(isnan(D10x_DE)))=0;
D10x_Norm(find(isnan(D10x_Norm)))=0;
mNR=mean(D60x_Norm);
vNR=std(D60x_Norm);
xa=[0.85 1.15];
bar(xa,[mNR,mDE])
errorbar([mNR,mDE],[vNR,vDE])
hold on
bar(xa,[mNR,mDE])
errorbar(xa,[mNR,mDE],[vNR,vDE])
hold on
bar(xa,[mNR,mDE])
hist(D60x_Norm)
hist(D60x_DE)
hist(D60x_Norm)
cumsum(D60x_DE)
h=hist(D60x_Norm,0:30:900)
bar(0:30:900,h)
h=hist(D60x_DE,0:30:900)
bar(0:30:900,h)
hcum=cumsum(h)
plot(0:30:900,hcum)
plot(0:30:900,hcum./33)
hold on
h=hist(D60x_Norm,0:30:900)
hcum=cumsum(h);
plot(0:30:900,hcum./31)
x=-10:1:10;
y=0:1:20;
randn(length(x),length(y))
M=randn(length(x),length(y))
find(M>0)
[x,y]=find(M>0)
x
y
doc arctan
test
doc arctan
arctan([1,2])
atan(0.5)
arctan(0.5)
atan(1)
doc atan
size(avg)
size(xi)
angy
angx
size(xi)
angy
angx
ScatterPlot_widthDistribution_darkexp
hist(Ang4_60,-3.14/2:0.1:3.14/2)
figure(1)
hist(Ang4_60,-3.14/2:0.1:3.14/2)
Ang4_60
ScatterPlot_widthDistribution_darkexp
Ang4_60
hist(Ang4_60,-3.14/2:0.1:3.14/2)
hist(Ang4_10,-3.14/2:0.1:3.14/2)
angx=-3.14/2:0.1:3.14/2;
x=sin(angx);
x=3sin(angx);
x=3.*sin(angx);
y=3.*cos(angx);
h=hist(Ang4_60,-3.14/2:0.1:3.14/2);
h=h./sum(h);
plot(h)
figure(2)
plot(h)
figure(2)
plot(h)
cc=jet(100);
hind=ceil(h.*100);
h=hist(Ang4_60,-3.14/2:0.1:3.14/2);
sum(h)
plot(h)
hist(Ang4_60,-3.14/2:0.1:3.14/2);
hist(Ang4_60,-3.14/2:0.2:3.14/2);
h=hist(Ang4_60,-3.14/2:0.1:3.14/2);
plot(x,y,'color',cc(h))
h=ceil(h)
plot(x,y,'color',cc(h))
plot(x,y,'color',cc(h+1))
plot(x,y,'color',cc(h)+1)
h(find(h==0))=1:
h(find(h==0))=1;
plot(x,y,'color',cc(h)+1)
plot(x,y,'color',cc(h))
h=ceil(h);
plot(x,y,'color',cc(:,h))
plot(x,y,'color',cc(h,:))
cc(1,:)
plot(x,y,'color',cc(h,:))
plot(x,y,'color',float(cc(h,:)))
plot(x,y,'color',double(cc(h,:)))
ScatterPlot_widthDistribution_darkexp
hist(Ang4_60,-3.14/2:0.2:3.14/2);
ScatterPlot_widthDistribution_darkexp
hist(Ang4_60,-3.14/2:0.2:3.14/2);
hist(Ang4_10,-3.14/2:0.2:3.14/2);
hist(Ang4_60,-3.14/2:0.2:3.14/2);
x=3.*sin(angx);
y=3.*cos(angx);
[X,Y]=meshgrid(x,y)
X
h=hist(Ang4_60,-3.14/2:0.1:3.14/2);
size(h)
size(X)
x
X
r=2:0.2:4;
X=r*x;
r=(2:1:length(x)+2);
size(r)
size(x)
r=(2:1:length(x)+1);
dr=2./length(x);
r=2:dr:4;
size(r)
r=r(1:end-1);
X=r*x;
size(x)
size(r)
X=r'*x;
size(X)
Y=r'*y;
size(Y)
pcolor(X,Y,h)
C=r'*h
pcolor(X,Y,C)
figure(100)
pcolor(X,Y,C)
C=ones(1,length(r))'*h
pcolor(X,Y,C)
h=hist(Ang4_60,-3.14/2:0.2:3.14/2);
angx=-3.14/2:0.2:3.14/2;
x=3.*sin(angx);
y=3.*cos(angx);
dr=2./length(x);
r=2+dr:dr:4;
r=2:1:4;
X=r'*x;
Y=r'*y;
C=ones(1,length(r))'*h
pcolor(X,Y,C)
h=hist(Ang4_60,-3.14/2:0.2:3.14/2);
h=h./sum(h);
C=ones(1,length(r))'*h
pcolor(X,Y,C)
h=hist(Ang4_10,-3.14/2:0.2:3.14/2);
C=ones(1,length(r))'*h
pcolor(X,Y,C)
h=h./sum(h);
C=ones(1,length(r))'*h
pcolor(X,Y,C)
ScatterPlot_widthDistribution_darkexp
h=hist(Ang4_60,-3.14/2:0.2:3.14/2);
h=h./sum(h);
C=ones(1,length(r))'*h
pcolor(X,Y,C)
size(X)
size(Y)
x=3.*sin(angx);
y=3.*cos(angx);
X=r'*x;
Y=r'*y;
size(X)
C=ones(1,length(r))'*h
pcolor(X,Y,C)
r=2:1:3;
X=r'*x;
Y=r'*y;
C=ones(1,length(r))'*h;
pcolor(X,Y,C)
ScatterPlot_widthDistribution_darkexp
h=hist(Ang4_60,-3.14/2:0.2:3.14/2);
h=h./sum(h);
X=r'*x;
Y=r'*y;
C=ones(1,length(r))'*h;
figure(2)
pcolor(X,Y,C)
figure(2)
pcolor(X,Y,C)
ScatterPlot_widthDistribution_darkexp
h=hist(Ang4_60,-3.14/2:0.2:3.14/2);
h=h./sum(h);
X=r'*x;
Y=r'*y;
C=ones(1,length(r))'*h;
figure(1)
pcolor(X,Y,C)
Ppk60=Ppeaklayer60;
Pchg60=Pchglayer60;
Ppk10=Ppeaklayer10;
Pchg10=Pchglayer10;
ScatterPlot_widthDistribution_darkexp
Ppk60DE=Ppeaklayer60;
Pchg60DE=Pchglayer60;
Ppk10DE=Ppeaklayer10;
Pchg10DE=Pchglayer10;
mPk60=mean(Ppk60);
vPk60=std(Ppk60);
mPk10=mean(Ppk10);
vPk10=std(Ppk10);
mPk60DE=mean(Ppk60DE);
vPk60DE=std(Ppk60DE);
mPk10DE=mean(Ppk10DE);
vPk10DE=std(Ppk10DE);
mchg60=mean(Pchg60);
vchg60=std(Pchg60);
mchg10=mean(Pchg10);
vchg10=std(Pchg10);
mchg60DE=mean(Pchg60DE);
vchg60DE=std(Pchg60DE);
mchg10DE=mean(Pchg10DE);
vchg10DE=std(Pchg10DE);
hist(mPk60)
figure(1)
hist(mPk60)
mPk60
mPk60DE
bar(mPk60DE)
ScatterPlot_widthDistribution_darkexp
mPk60=mean(Ppk60);
vPk60=std(Ppk60);
mPk10=mean(Ppk10);
vPk10=std(Ppk10);
mPk60
mPk60=nanmean(Ppk60);
vPk60=std(Ppk60);
mPk10=mean(Ppk10);
vPk10=std(Ppk10);
mPk60
mPk60=nanmean(Ppk60);
vPk60=nanstd(Ppk60);
mPk10=nanmean(Ppk10);
vPk10=nanstd(Ppk10);
Ppk60(find(isnan(Ppk60)))=0;
mPk60=mean(Ppk60);
vPk60=std(Ppk60);
mPk10=mean(Ppk10);
vPk10=std(Ppk10);
mPk60
ScatterPlot_widthDistribution_darkexp
mPk60=nanmean(Ppk60);
vPk60=nanstd(Ppk60);
mPk10=nanmean(Ppk10);
vPk10=nanstd(Ppk10);
Ppk60=Ppeaklayer60;
Pchg60=Pchglayer60;
Ppk10=Ppeaklayer10;
Pchg10=Pchglayer10;
mPk60=nanmean(Ppk60);
vPk60=nanstd(Ppk60);
mPk10=nanmean(Ppk10);
vPk10=nanstd(Ppk10);
mPk60
mchg60=nanmean(Pchg60);
vchg60=nanstd(Pchg60);
mchg10=nanmean(Pchg10);
vchg10=nanstd(Pchg10);
PK=[mPk60;mPk60DE];
Vpk=[vPk60;vPk60DE];
errorbar([1 2 3 4; 1.1 2.1 3.1 4.1],PK,Vpk)
errorbar([0.9 1.9 2.9 3.9; 1.1 2.1 3.1 4.1],PK,Vpk)
errorbar([0.8 1.8 2.8 3.8; 1.2 2.2 3.2 4.2],PK,Vpk)
hold on
bar([0.8 1.8 2.8 3.8; 1.2 2.2 3.2 4.2],PK)
figure(1)
bar([0.8 1.8 2.8 3.8; 1.2 2.2 3.2 4.2],PK)
figure(2)
errorbar([0.8 1.8 2.8 3.8; 1.2 2.2 3.2 4.2],PK,Vpk)
doc errorbar
bar([0.8 1.8 2.8 3.8; 1.2 2.2 3.2 4.2],PK)
hold on
x=[3.67 3.75 3.83 3.91 4.09 4.17 4.25 4.33];
y=[PK(1,:) PK(2,:)];
errorbar(x,y,Vpk(:))
x=[3.64 3.81 3.88 3.91 4.09 4.17 4.25 4.33];
errorbar(x,y,Vpk(:))
x=[3.73 3.79 3.85 3.91 4.09 4.17 4.25 4.33];
errorbar(x,y,Vpk(:))
x=[3.7 3.77 3.84 3.91 4.09 4.17 4.25 4.33];
errorbar(x,y,Vpk(:))
x=[3.685 3.76 3.835 3.91 4.09 4.17 4.25 4.33];
errorbar(x,y,Vpk(:))
x=[3.687 3.762 3.835 3.91 4.09 4.17 4.25 4.33];
errorbar(x,y,Vpk(:))
x=[3.688 3.762 3.835 3.91 4.09 4.17 4.25 4.33];
errorbar(x,y,Vpk(:))
x=[3.689 3.762 3.835 3.91 4.09 4.165 4.24 4.31];
errorbar(x,y,Vpk(:))
errorbar(x,y,Vpk(:),'k*')
bar([0.8 1.8 2.8 3.8; 1.2 2.2 3.2 4.2],PK)
y=[PK(1,:) PK(2,:)];
hold on
errorbar(x,y,Vpk(:),'k*')
figure(1)
errorbar(x,y,Vpk(:),'k*')
hold on
bar([0.8 1.8 2.8 3.8; 1.2 2.2 3.2 4.2],PK)
signrank
doc signrank
hist(Ppk60(:,3))
figure(100)
hist(Ppk60(:,3))
h=hist(Ppk60(:,3),[0:0.1:1]);
hist(Ppk60(:,3),[0:0.1:1]);
hDE=hist(Ppk60(:,3),[0:0.1:1]);
hDE=hist(Ppk60DE(:,3),[0:0.1:1]);
figure(101);hist(Ppk60DE(:,3),[0:0.1:1]);
[h,p]=lillietest(Ppk60)
[h,p]=lillietest(Ppk60(:,3))
doc lillietest
[h,p]=lillietest(Ppk60DE(:,3))
doc ranksum
[p,h]=ranksum(Ppk60(:,3),Ppk60DE(:,3))
[p,h]=ranksum(Ppk60(:,1),Ppk60DE(:,1))
[p,h]=ranksum(Ppk60(:,2),Ppk60DE(:,2))
[p,h]=ranksum(Ppk60(:,4),Ppk60DE(:,4))
bar([0.8 1.8 2.8 3.8; 1.2 2.2 3.2 4.2],PK)
hold on
errorbar(x,y,Vpk(:),'k*')
box off
[p4,h4]=ranksum(Ppk10(:,3),Ppk61DE(:,3))
[p1,h1]=ranksum(Ppk10(:,1),Ppk10DE(:,1))
[p2,h2]=ranksum(Ppk10(:,2),Ppk10DE(:,2))
[p5,h5]=ranksum(Ppk10(:,4),Ppk10DE(:,4))
[p4,h4]=ranksum(Ppk10(:,3),Ppk60DE(:,3))
[p1,h1]=ranksum(Ppk10(:,1),Ppk10DE(:,1))
[p2,h2]=ranksum(Ppk10(:,2),Ppk10DE(:,2))
[p5,h5]=ranksum(Ppk10(:,4),Ppk10DE(:,4))
PK=[mPk10;mPk10DE];
bar([0.8 1.8 2.8 3.8; 1.2 2.2 3.2 4.2],PK)
mPk10=nanmean(Ppk10);
mPk10DE=nanmean(Ppk10DE);
vPk10=nanstd(Ppk10);
vPk10DE=nanstd(Ppk10DE);
x=[3.689 3.762 3.835 3.91 4.09 4.165 4.24 4.31];
y=[PK(1,:) PK(2,:)];
Vpk=[vPk10 vPk10DE];
errorbar(x,y,Vpk(:),'k*')
figure(100)
errorbar(x,y,Vpk(:),'k*')
vPk10DE
vPk10
=nanmean(Ppk10DE);
mPk10DE
mPk10
PK=[mPk10;mPk10DE];
errorbar(x,y,Vpk(:),'k*')
x=[3.689 3.762 3.835 3.91 4.09 4.165 4.24 4.31];
y=[PK(1,:) PK(2,:)];
errorbar(x,y,Vpk(:),'k*')
figure(100)
errorbar(x,y,Vpk(:),'k*')
hold on
bar([0.8 1.8 2.8 3.8; 1.2 2.2 3.2 4.2],PK)
box off
RR0=RR;
gex_mean_Vwh(Vmean )
hold on
plot(RR0,Vmean)
gex_mean_Vwh(Vmean )
gex_mean_Vwh(Vmean,0)
hold on
plot(RR0,Vmean)
poissionrate_vmean_h
c
size(c)
image(c)
poissionrate_vmean_h
hold on
gex_mean_Vwh(Vmean,0)
poissionrate_vmean_h
gex_mean_Vwh(Vmean,10)
poissionrate_vmean_h
gex_mean_Vwh(Vmean,20)
firingrate_poisson_vwh_combine
poissionrate_vmean_h
hold on
gex_mean_Vwh
gex_mean_Vwh(Vmean,20 )
gex_mean_Vwh
gex_mean_Vwh(Vmean,0)
gex_mean_Vwh(Vmean,10)
gex_mean_Vwh(Vmean,20)
poissionrate_vmean_h
gex_mean_Vwh(Vmean,0)
gex_mean_Vwh
gex_mean_Vwh(Vmean,0)
gex_mean_Vwh
gex_mean_Vwh(Vmean,0)
poissionrate_vmean_h
gex_mean_Vwh(Vmean,10)
poissionrate_vmean_h
gex_mean_Vwh(Vmean,0)
poissionrate_vmean_h
gex_mean_Vwh(Vmean,20)
figure(100)
plot(RR,Fr)
firingrate_poisson_vw
%-- 3/5/14, 10:35 AM --%
firingrate_poisson_vw
%-- 3/10/14, 11:09 AM --%
test
%-- 3/11/14, 2:59 PM --%
test
1-NaN
nanmean(NaN)
test
mapave_gui2
a
Indx
a{376}
a(376)
a(372)
Ind
pt
C=load(fullfile(folderpath,a(bs).name));
C.cell_folde
C
C.cellls
cell_folder
cells.cell_folder
cells
C
C.cells
C.cells.cell_folder
C=load(fullfile(folderpath,a(Indx).name));
C.cells.cell_folder
C.cells.mapname
C.cells.cell_folder='/Users/lindameng/Documents/study/meng/ttx/111412/bs0004';
doc save
a(Indx)
a(Indx).name
filename=fullfile(folderpath,a(Indx).name)
cells=C.cells;
save(filename,'cells')
C=load(fullfile(folderpath,a(Indx).name));
C.cells
cells
load '/Volumes/disk1/study/meng_code2/silent/eventwindow50/Map_121114_1719'
load '/Volumes/disk1/study/meng_code2/silent/eventwindow50/LSPS_Map_121114_1719'
cells.cellfolder
cells
load '/Volumes/disk1/study/meng_code2/silent/eventwindow50/LSPS_Map_121114_1512'
cells
load '/Volumes/disk1/study/meng_code2/silent/eventwindow50/LSPS_Map_121114_1706'
cells
load '/Volumes/disk1/study/meng_code2/silent/eventwindow50/LSPS_Map_121114_1685'
cells
load '/Volumes/disk1/study/meng_code2/silent/eventwindow50/LSPS_Map_121114_1727'
cells
load '/Volumes/disk1/study/meng_code2/silent/eventwindow50/LSPS_Map_121114_1744'
cells
load '/Volumes/disk1/study/meng_code2/silent/eventwindow50/LSPS_Map_121114_1751'
load '/Volumes/disk1/study/meng_code2/silent/eventwindow50/LSPS_Map_121114_1852'
cells
load '/Volumes/disk1/study/meng_code2/silent/eventwindow50/LSPS_Map_121114_1901'
cells
load '/Volumes/disk1/study/meng_code2/silent/eventwindow50/LSPS_Map_121114_1920'
cells
cells.cell_folder='/Users/lindameng/Documents/study/meng/ttx/111412/bs0005'
save('/Volumes/disk1/study/meng_code2/silent/eventwindow50/LSPS_Map_121114_1920',cells)
filename='/Volumes/disk1/study/meng_code2/silent/eventwindow50/LSPS_Map_121114_1920';
save(filename,'cells')
load '/Volumes/disk1/study/meng_code2/silent/eventwindow50/LSPS_Map_121114_1929'
cells
mapave_gui2
a(378)
a(280)
a(380)
load '/Volumes/disk1/study/meng_code2/silent/eventwindow50/LSPS_Map_121114_2008'
cells.cell_folder='/Users/lindameng/Documents/study/meng/ttx/111412/bs0006'
filename='/Volumes/disk1/study/meng_code2/silent/eventwindow50/LSPS_Map_121114_2008';
save(filename,'cells')
mapave_gui2
a(384)
load '/Volumes/disk1/study/meng_code2/silent/eventwindow50/LSPS_Map_121126_1413'
cells
cells.cell_folder='/Users/lindameng/Documents/study/meng/ttx/112612/bs0001'
filename='/Volumes/disk1/study/meng_code2/silent/eventwindow50/LSPS_Map_121126_1413';
save(filename,'cells')
a(380)
a(382)
a(385)
a(384)
a(388)
a(389)
a(390)
a(388)
guide
ScatterPlot_widthDistribution
%-- 4/8/14, 9:25 AM --%
plot(RR,Vmean)
hold on
plot(RR,Vmean)
gex_mean_Vwh(Vmean,0)
gex_mean_Vwh(Vmean,10)
gex_mean_Vwh(Vmean,20)
box off
firingrate_poisson_vwh_combine
ScatterPlot_widthDistribution
test
ScatterPlot_widthDistribution
Cellrecord{1}
ScatterPlot_widthDistribution
Cellrecord{3}
ScatterPlot_widthDistribution
x=[0.5 0.4038,0.4048]
y=x;
x=[1 2 3]
plotregression(x,y)
doc plotregression
regstats(y,x)
stat=regstats(y,x)
stat.fstat
stat.R
stat.Q
stat.r
doc regstats
x=[x';ones(3,1)]
x=[1 1;2 1;3 1]
help regress
[B,BINT,R,RINT,STATS] = regress(y,x)
y
y=y';
[B,BINT,R,RINT,STATS] = regress(y,x)
y=[0.2692 0.25, 0.1190]
[B,BINT,R,RINT,STATS] = regress(y,x)
[B,BINT,R,RINT,STATS] = regress(y',x)
y=[0.1154;0.1346;0.1667]
[B,BINT,R,RINT,STATS] = regress(y',x)
[B,BINT,R,RINT,STATS] = regress(y,x)
plot(x(:,1),y)
figure(100)
plot(x(:,1),y)
y=[0.2692 0.25, 0.1190]
plot(x(:,1),y)
hold on
y=[0.1154;0.1346;0.1667]
[B,BINT,R,RINT,STATS] = regress(y,x)
plot(x(:,1),y)
y=[0.1154;0.2115;0.2857]
[B,BINT,R,RINT,STATS] = regress(y,x)
plot(x(:,1),y)
y=[0.5;0.4038;0.4048];
plot(x(:,1),y)
[B,BINT,R,RINT,STATS] = regress(y,x)
y=
y
figure(100)
plotregression(x(:,1),y)
hold on
plot(x(:,1),y,'b*-')
y=[0.2692 0.25, 0.1190]
plotregression(x(:,1),y)
hold on
plot(x(:,1),y,'b*-')
y=[0.5;0.4038;0.4048];
figure(100)
plotregression(x(:,1),y)
hold on
plot(x(:,1),y,'b*-')
[B,BINT,R,RINT,STATS] = regress(y,x)
plot(x(:,1),y,'b*-')
hold on
plot(x(:,1),B(1).*x(:,1)+B(2))
y=[0.2692 0.25, 0.1190]
[B,BINT,R,RINT,STATS] = regress(y,x)
[B,BINT,R,RINT,STATS] = regress(y',x)
plot(x(:,1),y,'b*-')
plot(x(:,1),B(1).*x(:,1)+B(2),'b')
y=[0.1154;0.2115;0.2857]
[B,BINT,R,RINT,STATS] = regress(y,x);
plot(x(:,1),y,'b*-')
plot(x(:,1),B(1).*x(:,1)+B(2),'b')
y=[0.1154;0.1346;0.1667]
[B,BINT,R,RINT,STATS] = regress(y,x);
plot(x(:,1),y,'b*-')
plot(x(:,1),B(1).*x(:,1)+B(2),'b')
%-- 4/23/14, 11:17 AM --%
test
ScatterPlot_widthDistribution
Tau50
ScatterPlot_widthDistribution
Tau60
Tau50
Ind60=find(Tau60>0)
mTau60=mean(Tau60(Ind60))
Ind50=find(Tau50>0)
mTau50=mean(Tau50(Ind50))
Ind60apv=find(Tau60apv>0);
Ind60apv=find(Tau60Apv>0);
mTau60apv=mean(Tau60Apv(Ind60apv));
Ind50apv=find(Tau50Apv>0);
mTau50apv=mean(Tau50Apv(Ind50apv));
mtau=[mTau60 mTau50 mTau60apv mTau50apv]./10;
vTau60apv=std(Tau60Apv(Ind60apv));
vTau50apv=std(Tau50Apv(Ind50apv));
vTau50=std(Tau50(Ind50));
vTau60=std(Tau60(Ind60));
vtau=[vTau60 vTau50 vTau60apv vTau50apv]./10;
errorbar([1 2 3 4],mtau,vtau)
Tau50apv
Tau50Apv
hist(Tau50apv)
hist(Tau50Apv)
hist(Tau50Apv(find(Tau50Apv>0)))
test
ScatterPlot_widthDistribution
Ind60=find(Tau60>0)
mTau60=mean(Tau60(Ind60))
Ind50=find(Tau50>0)
mTau50=mean(Tau50(Ind50))
Ind60apv=find(Tau60Apv>0);
mTau60apv=mean(Tau60Apv(Ind60apv));
Ind50apv=find(Tau50Apv>0);
mTau50apv=mean(Tau50Apv(Ind50apv));
vTau60apv=std(Tau60Apv(Ind60apv));
vTau50apv=std(Tau50Apv(Ind50apv));
vTau50=std(Tau50(Ind50));
vTau60=std(Tau60(Ind60));
vtau=[vTau60 vTau50 vTau60apv vTau50apv]./10;
mtau=[mTau60 mTau50 mTau60apv mTau50apv]./10;
errorbar([1 2 3 4],mtau,vtau)
Tau50apv
Tau50Apv
ScatterPlot_widthDistribution
Ind60=find(Tau60>0)
mTau60=mean(Tau60(Ind60))
Ind50=find(Tau50>0)
mTau50=mean(Tau50(Ind50))
Ind60apv=find(Tau60Apv>0);
mTau60apv=mean(Tau60Apv(Ind60apv));
Ind50apv=find(Tau50Apv>0);
mTau50apv=mean(Tau50Apv(Ind50apv));
vTau60apv=std(Tau60Apv(Ind60apv));
vTau50apv=std(Tau50Apv(Ind50apv));
vTau50=std(Tau50(Ind50));
vTau60=std(Tau60(Ind60));
vtau=[vTau60 vTau50 vTau60apv vTau50apv]./10;
mtau=[mTau60 mTau50 mTau60apv mTau50apv]./10;
errorbar([1 2 3 4],mtau,vtau)
Tau50apv
Tau50Apv
Tau50apv
Tau50Apv
hist(Tau50Apv)
hist(Tau50Apv(find(Tau50Apv<12000)))
Tau50Apv2=Tau50Apv(find(Tau50Apv<12000));
mTau50apv=mean(Tau50Apv2(Ind50apv));
Ind50apv=find(Tau50Apv2>0);
mTau50apv=mean(Tau50Apv2(Ind50apv));
mtau=[mTau60 mTau50 mTau60apv mTau50apv]./10;
vTau50apv=std(Tau50Apv2(Ind50apv));
vtau=[vTau60 vTau50 vTau60apv vTau50apv]./10;
errorbar([1 2 3 4],mtau,vtau)
hist(Tau60Apv)
b=deleteoutliers(Tau60Apv,0.05)
hist(b)
Ind60apv=find(b>0);
mTau60apv=mean(b(Ind60apv));
vTau60apv=std(b(Ind60apv));
b50=deleteoutliers(Tau50Apv)
Tau50Apv
mtau=[mTau60 mTau50 mTau60apv mTau50apv]./10;
vtau=[vTau60 vTau50 vTau60apv vTau50apv]./10;
errorbar([1 2 3 4],mtau,vtau)
hist(Tau60)
b50=deleteoutliers(Tau60)
Ind60=find(b>0);
mTau60=mean(b(Ind60));
vTau60=std(b(Ind60));
mtau=[mTau60 mTau50 mTau60apv mTau50apv]./10;
vtau=[vTau60 vTau50 vTau60apv vTau50apv]./10;
errorbar([1 2 3 4],mtau,vtau)
b
mean(b(find(b>0)))
mTau60
vtau=[vTau60 vTau50 vTau60apv vTau50apv]./10;
doc errorbar
vTau60
errorbar([1 2 3 4],mtau.*10,vtau.*10)
hist(Tau50)
b=deleteoutliers(Tau50,0.05)
hist(b)
b50=deleteoutliers(Tau60)
Ind50=find(Tau50>0)
mTau50=mean(b(Ind50))
Ind50=find(b>0)
mTau50=mean(b(Ind50))
vTau50=std(b(Ind50))
b=deleteoutliers(Tau60)
Ind60=find(b>0);
mTau60=mean(b(Ind60));
vTau60=std(b(Ind60));
b=deleteoutliers(Tau60Apv,0.05)
Ind60apv=find(b>0);
mTau60apv=mean(b(Ind60apv));
vTau60apv=std(b(Ind60apv));
b=deleteoutliers(Tau50Apv,0.05)
hist(Tau50Apv)
Tau50Apv
Ind50=find(Tau50Apv>0)
mtau=[mTau60 mTau50 mTau60apv mTau50apv];
vtau=[vTau60 vTau50 vTau60apv vTau50apv];
errorbar([1 2 3 4],mtau,vtau)
Tau50apv
Tau50Apv
Tau50Apv(findTau50Apv>0)
Tau50Apv(findTau50Apv>0))
Tau50Apv(find(Tau50Apv>0))
hist(Tau50Apv(find(Tau50Apv>0)))
t50apv=Tau50Apv(find(Tau50Apv>0));
b=deleteoutliers(t50Apv,0.05)
b=deleteoutliers(t50apv,0.05)
hist(b)
hist(b(b<10000))
mtau=[mTau60 mTau50 mTau60apv mean(b(b<10000))];
vtau=[vTau60 vTau50 vTau60apv std(b(b<10000))];
errorbar([1 2 3 4],mtau,vtau)
test
ScatterPlot_widthDistribution
Tau60
Tau60=Tau60(find(Tau60>0));
Tau50=Tau50(find(Tau50>0));
Tau60Apv=Tau60Apv(find(Tau60Apv>0));
Tau50Apv=Tau50Apv(find(Tau50Apv>0));
hist(Tau60)
figure(100)
hist(Tau60)
b60=deleteoutliers(Tau60,0.05);
b50=deleteoutliers(Tau50,0.05);
b60apv=deleteoutliers(Tau60Apv,0.05);
b50apv=deleteoutliers(Tau50Apv,0.05);
mTau60=mean(b60);
mTau50=mean(b50);
mTau60Apv=mean(b60apv);
mTau50Apv=mean(b50apv);
vTau60=std(b60);
vTau50=std(b50);
vTau60Apv=std(b60apv);
vTau50Apv=std(b50apv);
mtau=[mTau60 mTau50 mTau60Apv mTau50Apv];
vtau=[vTau60 vTau50 vTau60Apv vTau50Apv];
errorbar([1 2 3 4],mtau,vtau)
hist(b50apv)
ScatterPlot_widthDistribution
Tau60=Tau60(find(Tau60>0));
Tau50=Tau50(find(Tau50>0));
Tau60Apv=Tau60Apv(find(Tau60Apv>0));
Tau50Apv=Tau50Apv(find(Tau50Apv>0));
A60=A60(find(Tau60>0));
A50=A50(find(Tau50>0));
A60Apv=A60Apv(find(Tau60Apv>0));
A50Apv=A50Apv(Tau50Apv>0));
Tau60=Tau60(find(Tau60>0));
Tau50=Tau50(find(Tau50>0));
Tau60Apv=Tau60Apv(find(Tau60Apv>0));
Tau50Apv=Tau50Apv(find(Tau50Apv>0));
A60=A60(find(Tau60>0));
A50=A50(find(Tau50>0));
A60Apv=A60Apv(find(Tau60Apv>0));
A50Apv=A50Apv(Tau50Apv>0);
scatter(A60,Tau60)
scatter(log(A60),log(Tau60))
[b60,ind60]=deleteoutliers(Tau60,0.05);
[b50,ind50]=deleteoutliers(Tau50,0.05);
[b60apv,ind60a]=deleteoutliers(Tau60Apv,0.05);
[b50apv,ind50a]=deleteoutliers(Tau50Apv,0.05);
scatter(log(A60(b60)),log(b60))
scatter(log(A60(ind60)),log(b60))
ind60
b60
b60=deleteoutliers(Tau60,0.05);
b50=deleteoutliers(Tau50,0.05);
b60apv=deleteoutliers(Tau60Apv,0.05);
b50apv=deleteoutliers(Tau50Apv,0.05);
mTau60=mean(b60);
mTau50=mean(b50);
mTau60Apv=mean(b60apv);
mTau50Apv=mean(b50apv);
vTau60=std(b60);
vTau50=std(b50);
vTau60Apv=std(b60apv);
vTau50Apv=std(b50apv);
mtau=[mTau60 mTau50 mTau60Apv mTau50Apv];
vtau=[vTau60 vTau50 vTau60Apv vTau50Apv];
errorbar([1 2 3 4],mtau,vtau)
hist(b50apv)
mean(b50apv(b50apv<8000))
median(b50apv)
median(b50)
median(b60)
median(b60apv)
scatter(A60,Tau60)
scatter(A50Apv,Tau50Apv)
Napv
ScatterPlot_widthDistribution
hist(Tau60R)
Tau60R2=deleteoutliers(Tau60R,0.05);
hist(Tau60R2)
Tau50R2=deleteoutliers(Tau50R,0.05);
hist(Tau50R2)
mean(Tau50R2)
mean(Tau60R2)
std(Tau50R2)
std(Tau60R2)
Tau60R2
Tau50R2
size(Tau60R2)
ScatterPlot_widthDistribution
hist(SpInc)
hist(SpIncapv)
hist(SpIncApv)
mean(SpInc)
mean(SpIncApv)
size(SpInc)
b=deleteoutliers(SpIncApv,0.05)
mean(SpIncApv)
ScatterPlot_widthDistribution
Agapv
ScatterPlot_widthDistribution
Agapv
ScatterPlot_widthDistribution
(sum(Cellrecord{fn}.meanflagPtxInd50+Cellrecord{fn}.meanflagDPtxInd50)-sum(Cellrecord{fn}.meanflagPtxInd60+Cellrecord{fn}.meanflagDPtxInd60))
sum(Cellrecord{fn}.meanflagPtxInd60+Cellrecord{fn}.meanflagDPtxInd60)
ScatterPlot_widthDistribution
sum(Cellrecord{fn}.meanflagapvInd50+Cellrecord{fn}.meanflagDapvInd50)
sum(Cellrecord{fn}.meanflagapvInd60+Cellrecord{fn}.meanflagDapvInd60)
ScatterPlot_widthDistribution
SpIncApv
nanmean(SpIncApv)
nanmean(SpInc)
ScatterPlot_widthDistribution
Cellrecord{fn}.meanflagapvInd50+Cellrecord{fn}.meanflagDapvInd50)-sum(Cellrecord{fn}.meanflagapvInd60+Cellrecord{fn}.meanflagDapvInd60)
sum(Cellrecord{fn}.meanflagapvInd50+Cellrecord{fn}.meanflagDapvInd50)-sum(Cellrecord{fn}.meanflagapvInd60+Cellrecord{fn}.meanflagDapvInd60)
Cellrecord{fn}.meanflagapvInd50+Cellrecord{fn}.meanflagDapvInd50)
Cellrecord{fn}.meanflagapvInd50+Cellrecord{fn}.meanflagDapvInd50
Cellrecord{fn}.Pth
ScatterPlot_widthDistribution
SpInc
ScatterPlot_widthDistribution
SpInc
ScatterPlot_widthDistribution
SpInc
SpIncApv
SpInc1=SpInc(find(SpInc>0));
x1=ones(SpInc1,1);
x1=ones(length(SpInc1),1);
SpIncApv1=deleteoutliers(SpIncApv,0.05)
x2=ones(length(SpIncApv1),1).*2;
plot(x1,SpIncApv1,'o')
plot(x1,SpIncApv1','o')
plot(x1,SpInc,'o')
plot(x1,SpInc1,'o')
figure(100)
plot(x1,SpInc1,'o')
hold on
plot(x2,SpIncapv1,'o')
plot(x2,SpIncApv1,'o')
find(AgIn==8)
figure(101)
plot(ones(length(find(AgIn==8)),1),SpInc(find(AgIn==8)))
plot(ones(length(find(AgIn==8)),1),SpInc(find(AgIn==8)),'o')
hold on
plot(ones(length(find(Agapv==8)),1),SpIncapv(find(Agapv==8)),'o')
plot(ones(length(find(AgApv==8)),1),SpIncApv(find(AgApv==8)),'o')
plot(ones(length(find(AgApv==8)),1).*,SpIncApv(find(AgApv==8)),'o')
plot(ones(length(find(AgApv==8)),1).*2,SpIncApv(find(AgApv==8)),'o')
ScatterPlot_widthDistribution
Cellrecord{fn}
Cellrecord{fn}.meanflagDapvInd50
ScatterPlot_widthDistribution
SpIncApv
ScatterPlot_widthDistribution
SpIncApv
SpInc
cells.Tp
ScatterPlot_widthDistribution
SpInc
ScatterPlot_widthDistribution
Cellrecord{fn}
ScatterPlot_widthDistribution
Cellrecord{fn}
fn
ScatterPlot_widthDistribution
Cellrecord{1}
Cellrecord{2}
Cellrecord{3}
Cellrecord{4}
Cellrecord{5}
Cellrecord{6}
ScatterPlot_widthDistribution
x=ones(SpInc);
x=ones(size(SpInc));
x1=ones(size(SpInc));
x1=ones(size(SpInc)).*2;
x1=ones(size(SpInc));
x2=ones(size(SpIncApv));
plot(x1,SpInc,'ro')
hold on
plot(x2.*2,SpIncApv,'ko')
plot(x2.*2+randn(length(x2)).*0.1,SpIncApv,'ko')
doc randn
plot(x2.*2+randn(size(x2)).*0.1,SpIncApv,'ko')
plot(x2.*2+randn(size(x2)).*0.01,SpIncApv,'ko')
hold on
plot(x1,SpInc,'ro')
mP=mean(SpInc);
vP=std(SpInc);
mPapv=mean(SpIncApv);
vPapv=std(SpIncApv);
error([1,2],[mP,MPapv],[vP,vPapv])
errorbar([1,2],[mP,mPapv],[vP,vPapv])
errorbar([1,2],[mP,mPapv],[vP,vPapv],'.')
bar([1,2],[mP,mPapv])
box off
Cellrecord{fn}
ScatterPlot_widthDistribution
Nsi
minage=3;
maxage=5;
Ntotal=length(find(Ag>=minage&Ag<=maxage));
M35=[];
for i=1:Ntotal
end
Ind35=find(Ag>=minage&Ag<=maxage);
for i=1Ntotal
end
for i=1:Ntotal
M35=[M35;nansum(AreaCell60{i})];
M35_50=[M35_50;nansum(AreaCell50{i})]
M35si=[M35si;nansum(AreaCellsi{i})]
end
M35_50=[];
M35si=[];
M35=[];
for i=1:Ntotal
M35=[M35;nansum(AreaCell60{i})];
M35_50=[M35_50;nansum(AreaCell50{i})]
M35si=[M35si;nansum(AreaCellsi{i})]
end
doc accum
M69_50=[];
M69si=[];
M69=[];
Ntotal69=length(find(Ag>=6&Ag<=9));
Ind69=find(Ag>=minage&Ag<=maxage);
for i=1:Ntotal69
M69=[M69;nansum(AreaCell60{i})];
M69_50=[M69_50;nansum(AreaCell50{i})]
M69si=[M69si;nansum(AreaCellsi{i})]
end
M10_50=[];
M10si=[];
M10=[];
Ntotal10=length(find(Ag>=10&Ag<=15));
Ind10=find(Ag>=10&Ag<=15);
for i=1:Ntotal10
M10=[M10;nansum(AreaCell60{Ind10(i)})];
M10_50=[M10_50;nansum(AreaCell50{Ind10(i)})]
M10si=[M10si;nansum(AreaCellsi{Ind10(i)})]
end
M69_50=[];
M69si=[];
M69=[];
Ntotal69=length(find(Ag>=6&Ag<=9));
Ind69=find(Ag>=6&Ag<=9);
for i=1:Ntotal69
M69=[M69;nansum(AreaCell60{Ind69(i)})];
M69_50=[M69_50;nansum(AreaCell50{Ind69(i)})];
M69si=[M69si;nansum(AreaCellsi{Ind69(i)})];
end
M35_50=[];
M35si=[];
M35=[];
Ntotal35=length(find(Ag>=3&Ag<=5));
Ind35=find(Ag>=3&Ag<=5);
for i=1:Ntotal35
M35=[M35;nansum(AreaCell60{Ind35(i)})];
M35_50=[M35_50;nansum(AreaCell50{Ind35(i)})];
M35si=[M35si;nansum(AreaCellsi{Ind35(i)})];
end
cdfplot(M35)
hold on
cdfplot(M69)
hold on
cdfplot(M10)
box off
cdfplot(M35_50)
hold on
cdfplot(M69_50)
hold on
cdfplot(M10_50)
box off
cdfplot(M35si)
hold on
cdfplot(M69si)
hold on
cdfplot(M10si)
box off
[p,h]=ranksum(M35,M69)
[p,h]=ranksum(M35,M10)
[p,h]=ranksum(M69,M10)
[p,h]=ranksum(M35_50,M69_50)
[p,h]=ranksum(M35_50,M10_50)
[p,h]=ranksum(M69_50,M10_50)
[p,h]=ranksum(M35si,M69si)
[p,h]=ranksum(M35si,M10si)
[p,h]=ranksum(M69si,M10si)
%-- 5/14/14, 5:54 PM --%
firingrate_poisson_vwh_combine
gex_mean_Vwh(Vmean,0)
figure(1)
plot(Vmean,Fr)
hold on
gex_mean_Vwh(Vmean,10)
figure(1)
plot(Vmean,Fr,'r')
gex_mean_Vwh(Vmean,20)
figure(1)
plot(Vmean,Fr,'g')
box off
firingrate_poisson_vh_subtractive_varyInh
firingrate_poisson_vh_combine_varyInh
firingrate_poisson_vh_varyInh
clear
plot(RM FRM)
plot(RM,FRM)
i=1:length(RM);
i=1:length(RM)-1;
ind=find(FRM(i)<=100&FRM(i+1)>100)
ind=find(FRM(i)<=0.1&FRM(i+1)>0.1)
ind(3)
FRM(ind)
FRM(ind+1)
RM(ind)
ind=find(FRM(i)>=0.1&FRM(i+1)<0.1)
FRM(ind+1)
FRM(ind)
RM(ind+1)
ind=find(FRM(i)>=0.1&FRM(i+1)<0.1)
i=
i=1:length(RM)-1;
ind=find(FRM(i)>=0.1&FRM(i+1)<0.1)
RM(ind+1)
FRM(ind)
FRM(ind+1)
ind=find(FRM(i)<=0.1&FRM(i+1)>0.1)
FRM(ind+1)
FRM(ind)
RM(ind+1)
RM_vwh_Eular_STA_poison
RM_vwh_combine_STA_poison
InhEffect_vh
plot(RM(:,k),FRM(:,k))
InhEffect_vh
doc regress
InhEffect_vh
plot(Vmean,Fr)
plot(mean,Fr)
clear all
firingrate_poisson_vh_subtractive_varyInh
plot(RM(:,1),VmM(:,1))
plot(VmM(:,1),FRM(:,1))
hold on
plot(VmM(:,3),FRM(:,3))
plot(VmM(:,5),FRM(:,5))
box off
firingrate_poisson_vh_subtractive_varyEPSG
%-- 5/30/14, 12:36 PM --%
firingrate_poisson_vh_subtractive_varyEPSG
EGmax=0.5:0.5:10;
R=0.05:0.5:50;
[X,Y]=meshgrid(EGmax,R);
figure(1)
contourf(X,Y,FRM)
figure(2)
contourf(X,Y,VmM)
firingrate_poisson_vh_divisive_varyEPSG
firingrate_poisson_vh_subtractive_varyEPSG
EGmax=0.5:0.5:10;
R=0.05:0.5:50;
[X,Y]=meshgrid(EGmax,R);
figure(1)
contourf(X,Y,FRM)
figure(2)
contourf(X,Y,VmM)
%-- 6/10/14, 3:07 PM --%
firingrate_poisson_vh_combine_varyEPSG
firingrate_poisson_vh_subtractive_varyEPSG
firingrate_poisson_vh_combine_varyEPSG
EGmax=0.5:0.5:40;
R=0.05:0.5:50;
[X,Y]=meshgrid(EGmax,R);
figure(1)
contourf(X,Y,FRM)
figure(2)
contourf(X,Y,VmM)
firingrate_poisson_vh_combine_varyEPSG
%-- 9/9/14, 1:03 PM --%
firingrate_poisson_vh_subtractive_varyInh
%-- 9/18/14, 12:54 PM --%
doc parse_pv_pairs
doc get
selectdata
doc get
patch;surface;text;line
output = get(get(gca,'Children'),props)
get(gca,'Children')
axes
get(gca,'Children')
allchild(gca)
%-- 9/30/14, 2:28 PM --%
test
guide
%-- 10/6/14, 9:31 AM --%
test
doc boxplot
doc ecdf
load discrim
help strcat
load hospital
help strcat
help ordinal
C1 = [2 1 1 2;1 2 1 2;1 1 2 2;2 2 2 3]
C1 =
2     1     1     2
1     2     1     2
1     1     2     2
2     2     2     3
rank(C1)
T=cholcov(C1)
C2 = T'*T
C2.C1
help pdist
help gallery
C2 = T'*T
X = gallery('uniformdata',[10 2],12);
Y = pdist(X);
help pdist
Z = linkage(Y,'single');
clear
load fisheriris;
X = randn(150,10);
X(:,[1 3 5 7 ])= meas;
y = species;
c = cvpartition(y,'k',10);
help cvpartition
kstest
doc kstest
doc kstest2
doc ranksum
subplot
doc subplot
COMPARE_CONTROLVsObject
COMPARE_CONTROLVsObject('/Users/lindameng/Documents/study/kawens/Norm/meanValueeachcell_eventwindow50_direct8', '/Users/lindameng/Documents/study/kawens/DarkExp/meanValueeachcell_eventwindow50_direct8')
a
COMPARE_CONTROLVsObject('/Users/lindameng/Documents/study/kawens/Norm/meanValueeachcell_eventwindow50_direct8', '/Users/lindameng/Documents/study/kawens/DarkExp/meanValueeachcell_eventwindow50_direct8')
COMPARE_CONTROLVsObject('Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8','Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914')
a
folderMeanpath{cmpind}
a=dir(fullfile('Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8','*.mat'));
a
fullfile('Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8','*.mat')
doc dir
a=dir('Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8');
a
clear
y=dir
y=dir('/Volumes/disk1')
y=dir('/Volumes/disk1/Norm')
y=dir('/Volumes/disk1/kawens060914')
y=dir('/Volumes/disk1/kawens060914/Norm')
y=dir('/Volumes/disk1/kawens060914/Norm/meanValueeachcell')
y=dir('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8')
y=dir(fullfile('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8','8.mat'))
y=dir(fullfile('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8','*.mat'))
COMPARE_CONTROLVsObject('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914')
Cs
cs
Cs{1}
COMPARE_CONTROLVsObject
COMPARE_CONTROLVsObject('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914')
LdistPk60
COMPARE_CONTROLVsObject('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914')
Cg
Cg{1}
COMPARE_CONTROLVsObject('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914')
guide
clear
test
test
test
%-- 10/8/14, 3:07 PM --%
test
COMPARE_CONTROLVsObject('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914')
COMPARE_CONTROLVsObject
COMPARE_CONTROLVsObject('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914')
doc cdfplot
COMPARE_CONTROLVsObject
COMPARE_CONTROLVsObject('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914')
subplot(2,2,i)
cdfplot(Pchg601(:,i),'Color','r');
cdfplot(Pchg601(:,i),'c','r');
cdfplot(Pchg601(:,i),'r');
cdfplot(Pchg601(:,i));
set('Color',[0,0,0])
h=cdfplot(Pchg601(:,i))
set(h,'Color',[0,0,0])
COMPARE_CONTROLVsObject('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914')
COMPARE_CONTROLVsObject
COMPARE_CONTROLVsObject('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914')
COMPARE_CONTROLVsObject('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914')
doc ttest
COMPARE_CONTROLVsObject('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914')
COMPARE_CONTROLVsObject('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914',0)
COMPARE_CONTROLVsObject('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914',1)
COMPARE_CONTROLVsObject
COMPARE_CONTROLVsObject('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914',1)
close all
COMPARE_CONTROLVsObject('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914',1)
close all
COMPARE_CONTROLVsObject('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914',1)
close all
COMPARE_CONTROLVsObject('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914',1)
doc legend
COMPARE_CONTROLVsObject('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914',1)
close all
COMPARE_CONTROLVsObject('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914',1)
close all
COMPARE_CONTROLVsObject('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914',1)
close all
COMPARE_CONTROLVsObject('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914',1)
close all
COMPARE_CONTROLVsObject('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914',1)
doc xlim
COMPARE_CONTROLVsObject('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914',1)
COMPARE_CONTROLVsObject
COMPARE_CONTROLVsObject('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914',1)
COMPARE_CONTROLVsObject('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914',0)
%-- 10/10/14, 2:34 PM --%
firingrate_poisson_vh_subtractive_varyInh
figure
plot(T1,vv)
firingrate_poisson_vh_subtractive_varyInh
plot(T1,vv)
firingrate_poisson_vh_subtractive_varyInh
plot(T1,vv)
kh=6;%hinf(phasic kh=2; tonic kh=6)
delta_h=0;
v0=-63.619;
w0=(1+exp(-(v0+48)/6))^(-1/4);
m0=(1+exp(-(v0+38)/7))^(-1);
h0=(1+exp((v0+65+delta_h)/kh))^(-1);
h0
firingrate_poisson_vh_subtractive_varyInh
COMPARE_Maps
doc cells
COMPARE_Maps('/Users/xymeng/Documents/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8', '/Users/xymeng/Documents/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Users/xymeng/Documents/kawens060914')
%-- 10/13/14, 5:11 PM --%
COMPARE_Maps('/Users/xymeng/Documents/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8', '/Users/xymeng/Documents/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Users/xymeng/Documents/kawens060914')
%-- 10/13/14, 5:15 PM --%
COMPARE_Maps('/Users/xymeng/Documents/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8', '/Users/xymeng/Documents/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Users/xymeng/Documents/kawens060914')
COMPARE_Maps
COMPARE_Maps('/Users/xymeng/Documents/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8', '/Users/xymeng/Documents/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Users/xymeng/Documents/kawens060914')
COMPARE_Maps('/Users/xymeng/Documents/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8', '/Users/xymeng/Documents/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Users/xymeng/Documents/kawens060914')
%-- 10/13/14, 5:20 PM --%
COMPARE_Maps
COMPARE_Maps('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8', '/Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914')
%-- 10/13/14, 5:24 PM --%
COMPARE_Maps('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8', '/Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914')
COMPARE_Maps
COMPARE_Maps('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8', '/Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914')
%-- 10/15/14, 12:23 PM --%
COMPARE_Maps
COMPARE_Maps('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8', '/Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914')
doc interp2
COMPARE_Maps('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8', '/Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914')
doc interp2
COMPARE_Maps('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8', '/Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914')
close all
COMPARE_Maps('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8', '/Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914')
%-- 10/21/14, 9:22 AM --%
COMPARE_Maps('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8', '/Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914')
near neighbor
doc neighbor
doc gaps
gna_bar=1000;
kh=6;%hinf(phasic kh=2; tonic kh=6)
delta_h=6;
v0=-63.619;
w0=(1+exp(-(v0+48)/6))^(-1/4);
m0=(1+exp(-(v0+38)/7))^(-1);
h0=(1+exp((v0+65+delta_h)/kh))^(-1);
gna=360*0.11/h0
doc gscatter
doc classify
load fisheriris
ldaClass = classify(meas(:,1:2),meas(:,1:2),species);
[x,y] = meshgrid(4:.1:8,2:.1:4.5);
x = x(:);
y = y(:);
ScatterPlot_widthDistribution_HI
chg60N=chglayre60;
chg60N=chglayer60;
ScatterPlot_widthDistribution_HI
chg60N=chglayer60;
ScatterPlot_widthDistribution_HI
chg60D=chglayer60;
clusterN=2;
%     iii=find(~isnan(Pchg60)&~isnan(Pchg10)&~isnan(Pchgsi));
% Pchg60=Pchg60(iii);
% Pchg10=Pchg10(iii);
% Pchgsi=Pchgsi(iii);
% xext=xext(iii);
% yy=yy(iii);
X=[chg60N(:,3),chg60N(:,4)];
err=[];
Idcell=cell(1,100);
for ii=1:100
[Id,Centr,sumd]=kmeans(X,clusterN);
err=[err sum(sumd)];
Idcell{ii}=Id;
end
[Er iclust]=min(err);
Id=Idcell{iclust(1)};
%-- 10/22/14, 2:16 PM --%
ScatterPlot_widthDistribution_HI
chg60N=chglayer60;
chg60D=chglayer60;
ScatterPlot_widthDistribution_HI
chg60N=chglayer60;
doc classify
X=[chg60N(:,3),chg60N(:,4)];
err=[];
Idcell=cell(1,100);
for ii=1:100
[Id,Centr,sumd]=kmeans(X,clusterN);
err=[err sum(sumd)];
Idcell{ii}=Id;
end
[Er iclust]=min(err);
Id=Idcell{iclust(1)};
clusterN=2;
X=[chg60N(:,3),chg60N(:,4)];
err=[];
Idcell=cell(1,100);
for ii=1:100
[Id,Centr,sumd]=kmeans(X,clusterN);
err=[err sum(sumd)];
Idcell{ii}=Id;
end
[Er iclust]=min(err);
Id=Idcell{iclust(1)};
gscatter(X(find(Id==1),1)?X(find(Id==1),1))
gscatter(X(:,1),X(:,1),Id)
gscatter(X(:,1),X(:,1),Id)
gscatter(X(:,1),X(:,2),Id)
X2=[chg60D(:,3),chg60D(:,4)];
ldaClass = classify(X2,X,Id);
ldaClass
gscatter(X2(:,1),X2(:,2),ldaclass)
gscatter(X2(:,1),X2(:,2),ldaClass)
X=[chg60N(:,1),chg60N(:,2)?chg60N(:,3)];
err=[];
Idcell=cell(1,100);
for ii=1:100
[Id,Centr,sumd]=kmeans(X,clusterN);
err=[err sum(sumd)];
Idcell{ii}=Id;
end
[Er iclust]=min(err);
Id=Idcell{iclust(1)};
X=[chg60N(:,1),chg60N(:,2)?chg60N(:,3)];
err=[];
Idcell=cell(1,100);
for ii=1:100
[Id,Centr,sumd]=kmeans(X,clusterN);
err=[err sum(sumd)];
Idcell{ii}=Id;
end
[Er iclust]=min(err);
Id=Idcell{iclust(1)};X=[chg60N(:,1),chg60N(:,2),chg60N(:,3)];
err=[];
Idcell=cell(1,100);
for ii=1:100
[Id,Centr,sumd]=kmeans(X,clusterN);
err=[err sum(sumd)];
Idcell{ii}=Id;
end
[Er iclust]=min(err);
Id=Idcell{iclust(1)};
X=[chg60N(:,1),chg60N(:,2),chg60N(:,3)];
err=[];
Idcell=cell(1,100);
for ii=1:100
[Id,Centr,sumd]=kmeans(X,clusterN);
err=[err sum(sumd)];
Idcell{ii}=Id;
end
[Er iclust]=min(err);
Id=Idcell{iclust(1)};
scatter3(X(:,1),X(:,2),X(:,3))
scatter3(X(:,1),X(:,2),X(:,3),Id)
%-- 10/22/14, 4:51 PM --%
test
ScatterPlot_widthDistribution_HI
%-- 10/23/14, 6:52 PM --%
ScatterPlot_widthDistribution_HI
D60chgN=D60chg;
D10chgN=D10chgN;
D10chgN=D10chg;
ScatterPlot_widthDistribution_HI
D60chgD=D60chg;
D10chgD=D10chg;
clusterN=2;
X=[D60chgN(:,1),D60chgN(:,3)];
err=[];
Idcell=cell(1,100);
for ii=1:100
[Id,Centr,sumd]=kmeans(X,clusterN);
err=[err sum(sumd)];
Idcell{ii}=Id;
end
[Er iclust]=min(err);
Id=Idcell{iclust(1)};
Id
gscatter(X(:,1),X(:,2),Id)
X2=[D60chgD(:,1),D60chgD(:,3)];
ldaClass = classify(X2,X,Id);
gscatter(X2(:,1),X2(:,2),ldaclass)
gscatter(X2(:,1),X2(:,2),ldaClass)
gscatter(X(:,1),X(:,2),Id)
test
guide
X1=[X;X2]
err=[];
Idcell=cell(1,100);
for ii=1:100
[Id,Centr,sumd]=kmeans(X2,clusterN);
err=[err sum(sumd)];
Idcell{ii}=Id;
end
[Er iclust]=min(err);
Id=Idcell{iclust(1)};
err=[];
Idcell=cell(1,100);
for ii=1:100
[Id,Centr,sumd]=kmeans(X1,clusterN);
err=[err sum(sumd)];
Idcell{ii}=Id;
end
[Er iclust]=min(err);
Id=Idcell{iclust(1)};
gscatter(X1(:,1),X1(:,2),Id)
err=[];
Idcell=cell(1,100);
for ii=1:100
[Id,Centr,sumd]=kmeans(X1,3);
err=[err sum(sumd)];
Idcell{ii}=Id;
end
[Er iclust]=min(err);
Id=Idcell{iclust(1)};
gscatter(X1(:,1),X1(:,2),Id)
ldaClass = classify(X1,X1,Id);
gscatter(X2(:,1),X2(:,2),ldaclass)
gscatter(X2(:,1),X2(:,2),ldaClass)
gscatter(X2(:,1),X2(:,2),ldaClass)
gscatter(X1(:,1),X1(:,2),ldaClass)
x=0:1:1000;
y=0:1:1000;
XX=[x y];
size(XX)
XX=[x' y'];
size(XX)
gscldaClass = classify(XX,X1,Id);
gscatter(XX(:,1),XX(:,2),ldaClass)
gscatter(XX(:,1),XX(:,2),gscldaClass)
bad=-strcmp(ldaClass,Id)
[x,y]=meshgrid(x,y);
j= classify([x y],X1,Id);
size(X1)
size(Id)
x=x(:);
y=y(:);
j= classify([x y],X1,Id);
gscatter(x,y,j)
= classify([x y],meas(:,1:2),species)
doc fitcnb
help fitcnb
help statistictool
help statistic
hold on
gscatter(X1(:,1),X1(:,2),ldaClass)
load fisheriris
nbKD = fitcnb(meas(:,1:2), species, 'DistributionNames','kernel', 'Kernel','box');
nbKDResubErr = resubLoss(nbKD)
nbKDCV = crossval(nbKD, 'CVPartition',cp);
nbKDCVErr = kfoldLoss(nbKDCV)
labels = predict(nbKD, [x y]);
gscatter(x,y,labels,'rgb','osd')
D60chg=[D60chgN;D60chgD]?
gplotmatrix(D60chg,[],Cylinders,['c' 'b' 'm' 'g' 'r' 'y'],[],[],false)
gplotmatrix(D60chg,[],cylinder,['c' 'b' 'm' 'g' 'r' 'y'],[],[],false)
doc gplotmatrix
gplotmatrix(D60chg,[],cylinder)
size(X1)
size(D60chg)
size(ladClass)
size(ldaClass)
gplotmatrix(D60chg,[],ldaClass)
gscatter(x,y,j)
hold on
gscatter(X1(:,1),X1(:,2),ldaClass)
doc line
figure
line([1-0.2 1+0.2], [1 1])
err=[];
Idcell=cell(1,100);
for ii=1:100
[Id,Centr,sumd]=kmeans(X1,2);
err=[err sum(sumd)];
Idcell{ii}=Id;
end
[Er iclust]=min(err);
Id=Idcell{iclust(1)};
x=0:1:1000;
y=0:1:1000;
XX=[x y];
size(XX)
XX=[x' y'];
size(XX)
gscldaClass = classify(XX,X1,Id);
x=x(:);
y=y(:);
j= classify([x y],X1,Id);
gscatter(x,y,j)
[x,y]=meshgrid(x,y);
x=x(:);
y=y(:);
j= classify([x y],X1,Id);
gscatter(x,y,j)
hold on
gscatter(X1(:,1),X1(:,2),Id)
err=[];
Idcell=cell(1,100);
for ii=1:100
[Id,Centr,sumd]=kmeans(X1,2,'dist','sqeuclidean');
err=[err sum(sumd)];
Idcell{ii}=Id;
end
[silh2,h] = silhouette(X1,Id,'sqeuclidean');
gscatter(X1(:,1),X1(:,2),Id)
size(X1)
gscatter(X1(:,1),X1(:,2),Id)
err=[];
Idcell=cell(1,100);
for ii=1:100
[Id,Centr,sumd]=kmeans(X1,2,'dist','cos');
err=[err sum(sumd)];
Idcell{ii}=Id;
end
eucD=pdist(X1,'euclidean')
clustTreeEuc=likage(eucD,'average')
clustTreeEuc=linkage(eucD,'average')
cophenet(clustTreeEuc,eucD)
[h,nodes] = dendrogram(clustTreeEuc,0);
eucD
doc pdist
doc silhouette
Train_data = X1;
Result = [];
for num_of_cluster = 1:20
centroid = kmeans(Train_data,num_of_cluster,'distance','sqeuclid');
s = silhouette(Train_data,centroid,'sqeuclid');
Result = [ Result; num_of_cluster mean(s)];
end
plot( Result(:,1),Result(:,2),'r*-.');`
Train_data = X1;
Result = [];
for num_of_cluster = 1:10
centroid = kmeans(Train_data,num_of_cluster,'distance','sqeuclid');
s = silhouette(Train_data,centroid,'sqeuclid');
Result = [ Result; num_of_cluster mean(s)];
end
plot( Result(:,1),Result(:,2),'r*-.');`
Train_data = X1;
Result = [];
for num_of_cluster = 1:10
centroid = kmeans(Train_data,num_of_cluster,'distance','sqeuclid');
s = silhouette(Train_data,centroid,'sqeuclid');
Result = [ Result; num_of_cluster nanmean(s)];
end
plot( Result(:,1),Result(:,2),'r*-.');`
Train_data = X1;
Result = [];
for num_of_cluster = 1:10
centroid = kmeans(Train_data,num_of_cluster,'distance','sqeuclid');
s = silhouette(Train_data,centroid,'sqeuclid');
Result = [ Result; num_of_cluster nanmean(s)];
end
plot( Result(:,1),Result(:,2),'r*-.');
Train_data = X1;
Result = [];
for num_of_cluster = 1:10
centroid = kmeans(Train_data,num_of_cluster,'distance','correlation');
s = silhouette(Train_data,centroid,'correlation');
Result = [ Result; num_of_cluster nanmean(s)];
end
plot( Result(:,1),Result(:,2),'r*-.');
load ovariancancer;
whos
grp
ScatterPlot_widthDistribution_HI
chg60N=chglayer60;
chg10N=chglayer10;
ScatterPlot_widthDistribution_HI
chg60D=chglayer60;
chg10D=chglayer10;
doc factoran
%-- 10/30/14, 12:08 PM --%
x=[1,2,3];
y=[24 68 8];
plot(x,y,'*--r')
hold on
y=[48.1 36.5 15.4];
plot(x,y,'*--k')
box off
COMPARE_Maps('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8', '/Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914')
imagesc(pcaDens60)
figure
imagesc(pcaDens60)
COMPARE_Maps('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8', '/Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914')
figure
COMPARE_Maps('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8', '/Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914')
figure
imagesc(pcaDens60)
[xd,yd]=meshgrid(xi,yi)
xd
yd
dx=30;
indx=find(sqrt(xd.^2+yd.^2)<30-dx&&sqrt(xd.^2+yd.^2)>=30-dx);
indx=find(sqrt(xd.^2+yd.^2)<30-dx&sqrt(xd.^2+yd.^2)>=30-dx);
indx
indx=find(sqrt(xd.^2+yd.^2)<60-dx&sqrt(xd.^2+yd.^2)>=30);
indx
indx=find(sqrt(xd.^2+yd.^2)<30-dx&sqrt(xd.^2+yd.^2)>=60);
indx
xd
indx=find(sqrt(xd.^2+yd.^2)<60&sqrt(xd.^2+yd.^2)>=30);
indx
[indx indy]=find(sqrt(xd.^2+yd.^2)<60&sqrt(xd.^2+yd.^2)>=30);
indx
indy
avg60(indx,indy)
COMPARE_Maps('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8', '/Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914')
imagesc(pcaDens60)
figure()
imagesc(pcaDens60)
[indx indy]=find(sqrt(xd.^2+yd.^2)<60&sqrt(xd.^2+yd.^2)>=30);
avg60(indx indy)
avg60(indx,indy)
[indx indy]=find(sqrt(xd.^2+yd.^2)<90&sqrt(xd.^2+yd.^2)>=60);
avg60(indx,indy)
.2000
[indx indy]=find(sqrt(xd.^2+yd.^2)<120&sqrt(xd.^2+yd.^2)>=90);
avg60(indx,indy)
size(pcaDens60)
36*42
size(avg10(indx,indy))
[indx indy]
avg10(indx indy)
avg10(indx,indy)
avg10(indx(:),indy(:))
size(avg10(indx,indy))
for i=1:length(indx)
avg10(indx(i),indy(i))
end
j=1:length(indx);
avg10(indx(j),indx(j))
COMPARE_Maps('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8', '/Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914')
size(avg10(indx,indy))
size(pcaDens60)
36*42
COMPARE_Maps('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8', '/Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914')
size(pcaDens60)
COMPARE_Maps('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8', '/Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914')
size(pcaDens60)
36*42
figure()
imagesc(pcaDens60)
install_paths
test
doc cov
COMPARE_CONTROLVsObject
COMPARE_CONTROLVsObject('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8', '/Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914')
COMPARE_CONTROLVsObject('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8', '/Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914',1)
MeanValueCalculationSiRatio_HI
test
COMPARE_CONTROLVsObject('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8', '/Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914',1)
close all
COMPARE_CONTROLVsObject('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8', '/Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914',1)
COMPARE_CONTROLVsObject('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8', '/Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914',0)
close all
test
COMPARE_CONTROLVsObject('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8', '/Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914',0)
close all
COMPARE_CONTROLVsObject('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8', '/Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914',0)
doc boxplot
close all
COMPARE_CONTROLVsObject('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8', '/Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914',0)
doc sl2dpcaex
close all
COMPARE_CONTROLVsObject('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8', '/Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914',0)
COMPARE_Maps('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8', '/Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914')
COMPARE_Maps
COMPARE_Maps('/Volumes/disk1/kawens060914/Norm/meanValueeachcell_eventwindow50_direct8', '/Volumes/disk1/kawens060914/DarkExp/meanValueeachcell_eventwindow50_direct8','/Volumes/disk1/kawens060914')
[Mm, PL, PR] = sl2dpcaex(pcaDens60, [46,32], 102, 'dimfix')
size(pcaDens60)
[Mm, PL, PR] = sl2dpcaex(pcaDens60, [46;32], 102, 'dimfix')
[Mm, PL, PR] = sl2dpcaex(pcaDens60,[46,32], 102, 'dimfix')
[Mm, PL, PR] = sl2dpcaex(pcaDens60,[46,32], 102)
[Mm, PL, PR] = sl2dpcaex(pcaDens60,[46,32], 102,'grareduce')
doc sl2dpca_apply
Y = sl2dpca_apply(Mm, PL, PR, pcaDens60, [46,32], 102);
doc sl2dpca_construct
X = sl2dpca_construct(Mm, PL, PR, Y)
image(X(:,:,1))
figure
image(X(:,:,1))
image(X(:,:,2))
image(X(:,:,3))
image(X(:,:,4))
for i=1:102
imagesc(X(:,:,i))
end
close all
> for i=1:102
figure(i);imagesc(X(:,:,i))
end
for i=1:102
figure(i);imagesc(X(:,:,i))
end
size(Y)
for i=1:102
figure(i);imagesc(Y(:,:,i))
end
size(PL)
size(PR)
imagesc(X(:,:,i))
figure
imagesc(pcaDens60(:,:,i))
X1 = sl2dpca_construct(Mm, PL(:,1:2), PR(:,1:2), Y)
size(Y)
X1 = sl2dpca_construct(Mm, PL(:,1:2), PR(:,1:2), Y(1:2,1:2))
imagesc(X1(:,:,1))
X1 = sl2dpca_construct(Mm, PL(:,1:3), PR(:,1:3), Y(1:3,1:3))
imagesc(X1(:,:,1))
size(Y)
X1 = sl2dpca_construct(Mm, PL(:,1:3), PR(:,1:3), Y(1:3,1:3,:))
imagesc(X1(:,:,1))
X1 = sl2dpca_construct(Mm, PL(:,1:4), PR(:,1:4), Y(1:4,1:4,:))
imagesc(X1(:,:,1))
for i=1:102
subplot(2,1,1)
imagesc(pcaDens60(:,:,i))
subplot(2,1,2)
imagesc(X1(:,:,i))
end
close all
for i=1:102
subplot(2,1,1)
imagesc(pcaDens60(:,:,i))
subplot(2,1,2)
imagesc(X1(:,:,i))
end
for i=1:102
figure(i);subplot(2,1,1)
imagesc(pcaDens60(:,:,i))
subplot(2,1,2)
imagesc(X1(:,:,i))
end
Y(1:3,1:3,1)
PR(:,1:3)
PL(:,1:3)
doc trace
info
doc slsymeig
doc pca
Y(1:2,1:2)
size(Y)
Y(:,:,1)
X1 = sl2dpca_construct(Mm, PL(:,1), PR(:,1:3), Y(1,1:3,:))
for i=1:102
figure(i);subplot(2,1,1)
imagesc(pcaDens60(:,:,i))
subplot(2,1,2)
imagesc(X1(:,:,i))
end
figure
imagesc(Mm+PL(:,1)*(1,00)*PR(:,1:3))
size(PL)
imagesc(Mm+PL(:,1)*(1,00)*PR(:,1:3)')
imagesc(Mm+PL(:,1)*(1,0,0)*PR(:,1:3)')
imagesc(Mm+PL(:,1)*[1,0,0]*PR(:,1:3)')
imagesc(Mm+PL(:,1)*[0,1,0]*PR(:,1:3)')
imagesc(Mm+PL(:,1)*[0,0,1]*PR(:,1:3)')
imagesc(Mm+PL(:,2)*[0,0,1]*PR(:,1:3)')
imagesc(Mm+PL(:,3)*[0,0,1]*PR(:,1:3)')
imagesc(Mm+PL(:,4)*[0,0,1]*PR(:,1:3)')
imagesc(Mm+PL(:,5)*[0,0,1]*PR(:,1:3)')
imagesc(Mm+PL(:,6)*[0,0,1]*PR(:,1:3)')
imagesc(Mm+PL(:,1)*[0,0,1]*PR(:,1:3)')
imagesc(Mm+PL(:,1)*[0,1,0]*PR(:,1:3)')
imagesc(Mm+PL(:,1)*[1,0,0]*PR(:,1:3)')
imagesc(Mm+PL(:,1)*[0,1,0]*PR(:,1:3)')
imagesc(Mm+PL(:,1)*[1,0,0]*PR(:,1:3)')
X1 = sl2dpca_construct(Mm, PL(:,1), PR(:,1:3), Y(1,1:3,:))
X1 = sl2dpca_construct(Mm, PL(:,1:2), PR(:,1:3), Y(1:2,1:3,:))
for i=1:102
figure(i);subplot(2,1,1)
imagesc(pcaDens60(:,:,i))
subplot(2,1,2)
imagesc(X1(:,:,i))
end
X1 = sl2dpca_construct(Mm, PL(:,1:2), PR(:,1:2), Y(1:2,1:2,:))
for i=1:102
figure(i);subplot(2,1,1)
imagesc(pcaDens60(:,:,i))
subplot(2,1,2)
imagesc(X1(:,:,i))
end
X1 = sl2dpca_construct(Mm, PL(:,1:3), PR(:,1:3), Y(1:3,1:3,:))
for i=1:102
figure(i);subplot(2,1,1)
imagesc(pcaDens60(:,:,i))
subplot(2,1,2)
imagesc(X1(:,:,i))
end
X1 = sl2dpca_construct(Mm, PL(:,1), PR(:,1:6), Y(1,1:6,:))
for i=1:102
figure(i);subplot(2,1,1)
imagesc(pcaDens60(:,:,i))
subplot(2,1,2)
imagesc(X1(:,:,i))
end
X1 = sl2dpca_construct(Mm, PL(:,:), PR(:,1:3), Y(:,1:3,:))
for i=1:102
figure(i);subplot(2,1,1)
imagesc(pcaDens60(:,:,i))
subplot(2,1,2)
imagesc(X1(:,:,i))
end