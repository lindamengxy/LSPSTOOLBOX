  
folderMeanpath='/Users/lindameng/Documents/Dan/meanValueeachcell_eventwindow50_direct8';
a=dir(fullfile(folderMeanpath,'*.mat'));
Cellfile=size(a);

ZPeak=[];
ZChg=[];
ZDen=[];
Nsi=0;
X=[];
Y=[];
for i= 1:1:Cellfile
    Cg= load(fullfile(folderMeanpath,a(i).name));
    if isfield(Cg,'Cellrecord')
        Cellrecord=Cg.Cellrecord;
        for fn=1:Cg.Num_cell
             if isfield(Cellrecord{fn},'meanpeakttxInd60')
                 
                 Nsi=Nsi+1;
                 Soma= Cellrecord{fn}.SomaCoordinates;
                 flipimg=Cellrecord{fn}.flipimg;
                 pth=char(Cellrecord{fn}.Pth);
                 stimcoordinates=Cellrecord{fn}.StimCoordinates;
                 rotate_angle=Cellrecord{fn}.SpatialRotation;
                 [stimCoordinates]=Matrixrotate_singlemap_GUI(Soma,stimcoordinates,flipimg,rotate_angle);
                 X=[X, stimCoordinates(1,:)];
                 Y=[Y, stimCoordinates(2,:)];
%                  Xaxis{Nsi}=X;
%                  Yaxis{Nsi}=Y;
                 Zpk60=zeros(size(stimCoordinates(1,:))).*NaN;
                 Zchg60=zeros(size(stimCoordinates(1,:))).*NaN;
                 Zden60=zeros(size(stimCoordinates(1,:)));
                 
                 
                 Indden=find((Cellrecord{fn}.meanflagttxInd60>0));
                 Zden60(Indden)=Cellrecord{fn}.meanflagttxInd60(Indden);
                
                 Zpk60(Indden)=Cellrecord{fn}.meanpeakttxInd60(Indden);
                
                 Zchg60(Indden)=Cellrecord{fn}.meanareattxInd60(Indden);
                 
                 ZPeak=[ZPeak Zpk60];
ZChg=[ZChg Zchg60];
ZDen=[ZDen Zden60];

                 
             end
                 
        end
    end
end
                    
     Ncell=Nsi;  
            
 density=0.15;           
dens=min(10,0.8*Ncell);
    Xavg=X;
    Yavg=Y;
    dx=40;
    dy=40;
      xlim = [-500-dx-0.1 max(X)+dx];
    ylim = [-500-dy-0.1 500+dy];
    
    xtick = xlim(1):dx:xlim(2);
    ytick = ylim(1):dy:ylim(2);
    xi = xtick(1:end-1)+dx/2;
    yi = ytick(1:end-1)+dy/2;
    xii=xi(1):1:xi(end);
    yii=yi(1):1:yi(end);   
    
    [cnt_exc_dir ind] =griddensity2(Xavg,Yavg,[],[],xtick,ytick,Ncell,dens);
    
    
    sel = (cnt_exc_dir<=dens);
    isel = find(sel);
    selrm = [];
    for i=1:length(isel)
        sel2 = find(ind==isel(i)); selrm = [selrm sel2];
    end
    Xavg(selrm) = [];
     Yavg(selrm) = [];
    ZPeak(selrm) = [];
    ZChg(selrm) = [];
    ZDen(selrm) =[];
    
    [yi1,xi1]=meshgrid(yi,xi);
    [yii,xii]=meshgrid(yii,xii);
    avg = gridavg2(Xavg,Yavg,ZDen,xtick,ytick,false,Ncell,0);
    Ind2=find(avg<=density); 
    
    currfig= find_figure('Avgdensity60');
    clf;
    avg=avg';
    
    avg2=interp2(yi1,xi1,avg,yii,xii,'cubic');
    axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
        'DataAspectRatio',[1 1 1]);
    imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
    
    colorbar
    
    title(sprintf('#cell %i',Ncell))
    
    % charge
    avg = gridavg2(Xavg,Yavg,ZChg,xtick,ytick,false,Ncell,0);
    avg(Ind2)=0;
    
    currfig=  find_figure('Avgcharge60');
    clf;
    avg=avg';
    
    avg2=interp2(yi1,xi1,avg,yii,xii,'cubic');
    axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
        'DataAspectRatio',[1 1 1]);
    imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
    
    colorbar
    
    title(sprintf('#cell %i',Ncell))
    
    % Peak
       
    % charge
    avg = gridavg2(Xavg,Yavg,ZPeak,xtick,ytick,false,Ncell,0);
     avg(Ind2)=0;
    
    currfig=  find_figure('AvgPeak60');
    clf;
    avg=avg';
    
    avg2=interp2(yi1,xi1,avg,yii,xii,'cubic');
    axes1 = axes('Parent',currfig,'YDir','reverse','Layer','top',...
        'DataAspectRatio',[1 1 1]);
    imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;
   
    colorbar
    
    title(sprintf('#cell %i',Ncell))
   