function handles= MeanValueCalculationSiRatio( handles )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
datafolder=handles.avgfolder;
if isfield(handles,'folderMeanpath')
    folderMeanpath=handles.folderMeanpath;
else
    folderMeanname=sprintf('meanValueeachcell_eventwindow%i_direct%i',str2double(handles.eventwindow.String),str2double(handles.direct_t.String));
    folderMeanpath=fullfile(datafolder,folderMeanname);
end
a=dir(fullfile(folderMeanpath,'*.mat'));

Cellfile=size(a);
if ~isfield(handles,'avgTp')||isempty(handles.avgTp)
    handles.avgTp=1;
    set(handles.error_run,'string','Forget to choose the drugs type. The default version is High Mg. If not, please choose the right one')
else
    set(handles.error_run,'string','')
end
if ~isfield(handles,'minageValue')||isempty(handles.minageValue)
    handles.minage=0;
    
end
if ~isfield(handles,'cutoffdensityValue')||isempty(handles.cutoffdensityValue)
    handles.cutoffdensity=0;
    
end
if ~isfield(handles,'maxageValue')||isempty(handles.maxageValue)
    handles.maxage=Inf;
end
if ~isfield(handles,'maxdistValue')||isempty(handles.maxdistValue)
    handles.maxdist=Inf;
end
if ~isfield(handles,'mindistValue')||isempty(handles.mindistValue)
    handles.mindist=0;
    
end

if ~isfield(handles,'avgdarkexpValue')||isempty(handles.avgdarkexpValue)
    handles.avgdarkexp=0;
    
end
if ~isfield(handles,'avgHoldingPotentialValue')||isempty(handles.avgHoldingPotentialValue)
    handles.avgHoldingPotential=-60;
    set(handles.error_run,'string','Forget to fill in the holding potential. The default value is -60. If not, please choose the right one')
else
    set(handles.error_run,'string','')
end
if ~isfield(handles,'indb')||isempty(handles.indb)
   
    set(handles.error_run,'string','Forget to pick up the layer where the patched cell is')
else
    set(handles.error_run,'string','')
end


xaxis=[];
yaxis=[];
Pk=[];
Chg=[];
Den=[];
Lat=[];
dX=[];
AGE=[];
BD=[];
N=0;
Pchglayer=[];
Ppeaklayer=[];
Parealayer=[];
DistWidthPeak=[];
MCHG=[];
MPK=[];
DistWidthCharge=[];
SkPeaklayer=[];
Skchglayer=[];
Skdenslayer=[];
DistWidthDense=[];
DirArea=[];
indb=handles.indb;

%silent
xaxissi=[];
yaxissi=[];
Pksi=[];
Chgsi=[];
Densi=[];
Latsi=[];
dXsi=[];
AGEsi=[];
BDsi=[];
Nsi=0;
Pchgsilayer=[];
Ppeaksilayer=[];
Pareasilayer=[];
DistWidthPeaksi=[];
DistWidthChargesi=[];
AR=[];
DistWidthdenssi=[];
Perc=0.20;

for i= 1:1:Cellfile
     if ~strncmp(a(i).name,'.',1)
         a(i).name
    Cg= load(fullfile(folderMeanpath,a(i).name));
    
    if isfield(Cg,'Cellrecord')
    Cellrecord=Cg.Cellrecord;
    for fn=1:Cg.Num_cell
        
        h=Cellrecord{fn}.experimenttype;
         Bnd=Cellrecord{fn}.Boundry;
         sb=length(Bnd);
        % if handles.indb~=1
         if handles.indb<2
         Bnd1(1)=Bnd(1);
         Bnd1(2:4)=Bnd(3:5);  
         Bnd1(1,5)=Bnd(1,sb);
         elseif sb>6
             Bnd1(1)=Bnd(1);
         Bnd1(2:4)=Bnd(2:4);
         Bnd1(1,5)=Bnd(1,sb);
         elseif handles.indb==2
              Bnd1(1)=Bnd(1);
         Bnd1(2:4)=Bnd(3:5);  
         
             else Bnd1=Bnd;
         end
         Bnd=Bnd1;
         %end
         
%                   if size(Bnd,2)>3
%                       Bnd=Bnd(:,1:3);
%                   end
% %          Bnd1=zeros(1,size(Bnd,2)+1);
% %          Bnd1(1:4)=Bnd(1:4);
% %          Bnd1(6:end)=Bnd(5:end);
% %          Bnd1(5)=(Bnd(4)+Bnd(5))./2;
%           Bnd=Bnd1;
         AGE=[AGE,Cellrecord{fn}.age];
                BD=[BD;Bnd];
%          if isfield(Cellrecord{fn},'meanflagPtxInd60')&&isfield(Cellrecord{fn},'meanflagPtxInd50')
%              mflagPtxInd50=Cellrecord{fn}.meanflagPtxInd50;
%              mflagPtxInd60=Cellrecord{fn}.meanflagPtxInd60;
%              Matrix_NMDAonly=mflagPtxInd50.*~Cellrecord{fn}.meanflagPtxInd60.*~Cellrecord{fn}.meanflagDPtxInd60;
%              mratePtxInd50=Cellrecord{fn}.meanratePtxInd50;
%              Matrix_NMDAonly23=(mflagPtxInd50.*(mratePtxInd50>2/3)).*~(Cellrecord{fn}.meanflagPtxInd60.*(Cellrecord{fn}.meanratePtxInd60>2/3)).*~Cellrecord{fn}.meanflagDPtxInd60;
%              if nansum(mflagPtxInd50)
%                  AreaRatio=nansum(mflagPtxInd60)/nansum(mflagPtxInd50);
%              else
%                  AreaRatio=NaN;
%              end
%              AR=[AR,AreaRatio];
%          else 
%              Matrix_NMDAonly=[];
%              Matrix_NMDAonly23=[];
%              AreaRatio=NaN;
%          end
%          Cellrecord{fn}.Matrix_NMDAonly=Matrix_NMDAonly;
%         Cellrecord{fn}.Matrix_NMDAonly23=Matrix_NMDAonly23;
%         Cellrecord{fn}.AreaRatio= AreaRatio;
        if Cellrecord{fn}.age>=handles.minage&&Cellrecord{fn}.age<=handles.maxage&&abs(Cellrecord{fn}.darkexp-handles.avgdarkexp)<Perc&&strcmp(h,'Whole cell')
            
       
       
                
            if handles.avgTp==1&&abs(handles.avgHoldingPotential+60)<1&&isfield(Cellrecord{fn},'meanlatencyHighMg60')&&~isempty(Cellrecord{fn}.meanlatencyHighMg60)
                
                Soma= Cellrecord{fn}.SomaCoordinates;
                flipimg=Cellrecord{fn}.flipimg;
                pth=char(Cellrecord{fn}.Pth)
                stimcoordinates=Cellrecord{fn}.StimCoordinates;
                rotate_angle=Cellrecord{fn}.SpatialRotation;
                [stimCoordinates]=Matrixrotate_singlemap_GUI(Soma,stimcoordinates,flipimg,rotate_angle);
                X=stimCoordinates(1,:);
                Y=stimCoordinates(2,:);
                
                Zpk=zeros(size(X)).*NaN;
                Zchg=zeros(size(X)).*NaN;
                Zlat=zeros(size(X)).*NaN;
                Ymat=zeros(size(X)).*NaN;
                Xmat=zeros(size(X)).*NaN;
                if abs(Bnd(indb)-Bnd(1))>=handles.mindist&&abs(Bnd(indb)-Bnd(1))<=handles.maxdist
                    N=N+1
      
                    xaxis=[xaxis,X];
                    yaxis=[yaxis,Y];
                    dX=[dX,Cellrecord{fn}.PatternSpacing(1)];  
                    Zflag=Cellrecord{fn}.meanflagHighMg60;
                    
                    IndZ=find(Zflag>0); 
                    dist=sqrt(X(IndZ).^2+Y(IndZ).^2);
                    dist=sort(dist);
                    if length(dist)
                    Dist8060=dist(ceil(length(dist)*0.8));
                    else Dist8060=0;
                    end
                    Zpk(IndZ)=Cellrecord{fn}.meanpeakHighMg60(IndZ);          
                    Zchg(IndZ)=Cellrecord{fn}.meanareaHighMg60(IndZ);                 
                    Zlat(IndZ)=Cellrecord{fn}.meanlatencyHighMg60(IndZ);
                    Ymat(IndZ)=Y(IndZ);
                         Xmat(IndZ)=X(IndZ);
                         
                          Xmat=Xmat./sqrt(Xmat.^2+Ymat.^2);
                                Ymat=Ymat./sqrt(Xmat.^2+Ymat.^2);
                    ZflagD=Cellrecord{fn}.meanflagDHighMg60;
                    IndZD=find(ZflagD>0); 
                    Cellrecord{fn}.DirectAreaTp1Hp60=length(IndZD).*Cellrecord{fn}.PatternSpacing(1).^2;
                    DirArea=[DirArea;length(IndZD).*Cellrecord{fn}.PatternSpacing(1).^2];
                    
                    Pchg=[];
                    Ppeak=[];
                    Parea=[];
                    chg=[];
                    peak=[];
                    Mchg=[];
                    Mpeak=[];
                     Medianchg=[];
                    Medianpeak=[];
                    A60=[];
                    AxisSum60=[];
                    expSum60=[];
                    ChgexpSum60=[];
                    PkexpSum60=[];
                    Dist8060Layer=[];
                    for bi=1:length(Bnd)-1
                             if bi==1
                                 indBnd=find(X<Bnd(bi+1));
                             elseif bi==length(Bnd)-1
                                 indBnd=find(X>=Bnd(bi));
                             else
                            indBnd=find(X>=Bnd(bi)&X<Bnd(bi+1));
                        end
                        
                        Pchg=[Pchg nansum(Zchg(indBnd))/nansum(Zchg)];
                        Ppeak=[Ppeak nansum(Zpk(indBnd))/nansum(Zpk)];
                        Parea=[Parea nansum(Zflag(indBnd))/nansum(Zflag)];
                        chg=[chg nansum(Zchg(indBnd))];
                        peak=[peak nansum(Zpk(indBnd))];
                        Mchg=[Mchg nanmean(Zchg(indBnd))];
                        Mpeak=[Mpeak nanmean(Zpk(indBnd))];
                        Medianchg=[Medianchg nanmedian(Zchg(indBnd))];
                        Medianpeak=[Medianpeak nanmedian(Zpk(indBnd))];
                        A60=[A60 nansum(Zpk(indBnd)>0).*Cellrecord{fn}.PatternSpacing(1).^2];
                        AxisSum60=[AxisSum60 nansum(Ymat(indBnd))];
                        expSum60=[expSum60 abs(nansum((Xmat(indBnd)+i.*Ymat(indBnd))./sqrt(Xmat(indBnd).^2+Ymat(indBnd).^2)))./nansum(abs((Xmat(indBnd)+i.*Ymat(indBnd))./sqrt(Xmat(indBnd).^2+Ymat(indBnd).^2)))];
                        ChgexpSum60=[ChgexpSum60 abs(nansum(Zchg(indBnd).*(Xmat(indBnd)+i.*Ymat(indBnd))./sqrt(Xmat(indBnd).^2+Ymat(indBnd).^2)))./nansum(abs(Zchg(indBnd).*(Xmat(indBnd)+i.*Ymat(indBnd))./sqrt(Xmat(indBnd).^2+Ymat(indBnd).^2)))];
                        PkexpSum60=[PkexpSum60 abs(nansum(Zpk(indBnd).*(Xmat(indBnd)+i.*Ymat(indBnd))./sqrt(Xmat(indBnd).^2+Ymat(indBnd).^2)))./nansum(abs(Zpk(indBnd).*(Xmat(indBnd)+i.*Ymat(indBnd))./sqrt(Xmat(indBnd).^2+Ymat(indBnd).^2)))];
                        
                        indB=find(X(IndZ)>=Bnd(bi)&X(IndZ)<Bnd(bi+1));
                        if length(indB)>0
%                         dist=sqrt(X(IndZ(indB)).^2+Y(IndZ(indB)).^2);
                         dist=abs(Y(IndZ(indB)));
                        dist=sort(dist);
                         %Dist8060Layer=[Dist8060Layer abs(dist(ceil(length(dist)*0.1))-dist(ceil(length(dist)*0.9)))];
                        Dist8060Layer=[Dist8060Layer dist(ceil(length(dist)*0.8))];
                        else
                            Dist8060Layer=[Dist8060Layer 0];
                        end
                    end
                    
                    Pchglayer=[Pchglayer;Pchg];
                    Ppeaklayer=[Ppeaklayer;Ppeak];
                    Parealayer=[Parealayer;Parea];
                    
                    Pk=[Pk,Zpk];
                    Chg=[Chg, Zchg];
                    Lat=[Lat,Zlat];
                    Den=[Den Zflag'];
                    MCHG=[MCHG Mchg];
                    MPK=[MPK Mpeak];
                    
                    Cellrecord{fn}.PchgTp1Hp60=Pchg;
                    Cellrecord{fn}.PpeakTp1Hp60=Ppeak;
                    Cellrecord{fn}.PareaTp1Hp60=Parea;
                    Cellrecord{fn}.chgTp1Hp60=chg;
                    Cellrecord{fn}.peakTp1Hp60=peak;
                    Cellrecord{fn}.AreaTp1Hp60=A60;
                     Cellrecord{fn}.MchgTp1Hp60=Mchg;
                    Cellrecord{fn}.MpeakTp1Hp60=Mpeak;
                    Cellrecord{fn}.MedianchgTp1Hp60=Medianchg;
                    Cellrecord{fn}.MedianpeakTp1Hp60=Medianpeak;
                    Cellrecord{fn}.AxisSumTp1Hp60=AxisSum60;
                    Cellrecord{fn}.expSum60=expSum60;
                    Cellrecord{fn}.ChgexpSum60=ChgexpSum60;
                    Cellrecord{fn}.PkexpSum60=PkexpSum60;
                    Cellrecord{fn}.Dist8060Tp1Hp60=Dist8060;
                     Cellrecord{fn}.Dist8060LayerTp1Hp60=Dist8060Layer;
                    
                    % Peak---------------------------------------------------------------------------------------------------
                    dx =Cellrecord{fn}.PatternSpacing(1); dy = dx;
                    xlim = [ceil(min(X))-dx-0.1 ceil(max(X))]; ylim = [ceil(min(Y))-dy-0.1 ceil(max(Y))];
                    
                    
                    xtick = xlim(1):dx:xlim(2); ytick = ylim(1):dy:ylim(2);
                    xi = xtick(1:end-1)+dx/2; yi = ytick(1:end-1)+dy/2;
                    avg = gridavg2(X,Y,Zpk,xtick,ytick,false,1,0);
                    xx=nansum(avg);
                    yy=nansum(avg,2);
                    [mxx ii]=max(xx);
                    jj=find(xx>mxx*Perc);
                    
                    if ~isempty(jj)
                        %                      ctext=abs(xi(jj(1))-xi(jj(end)));
                        ctextratio=abs(xi(jj(1))-xi(jj(end)))/abs(max(X)-min(X));
                        Dx=max(dx,abs(xi(jj(1))-xi(jj(end))));
                    else Dx=0;ctext=0;ctextratio=0;
                    end
                    [myy ii]=max(yy);
                    jj=find(yy>myy*Perc);
                    if ~isempty(jj)
                        Dy=max(dx,abs(yi(jj(1))-yi(jj(end))));
                    else Dy=0;
                    end
                  
                    D=[];Skpeak=[];Laydist=[];
                    for bi=1:length(Bnd)-1
                        if bi==length(Bnd)-1
                             indBnd=find(xi>=Bnd(bi));
                        else
                            indBnd=find(xi>=Bnd(bi)&xi<Bnd(bi+1));
                        end
                        
                        y1=nansum(avg(:,indBnd),2);
                        s1=skewness(y1);
                        Laydist=[Laydist y1];
                     
                        [mx ii]=max(y1);
                        jj=find(y1>mx*Perc);
                        
                        if ~isempty(jj)
                            D1=max(dx,abs(yi(jj(1))-yi(jj(end))));
                        else D1=0;
                        end
                        D=[D D1];
                        Skpeak=[Skpeak s1];
                    end
                   
                    DW=[Dx ctextratio  Dy Dy./1000 D];
                    DistWidthPeak=[DistWidthPeak;DW];
                    SkPeaklayer=[SkPeaklayer;Skpeak];
                    Cellrecord{fn}.LaydistPkTp1Hp60=Laydist;
                    Cellrecord{fn}.DistWidthPeakTp1Hp60=DW;
                    Cellrecord{fn}.SkPeaklayerTp1Hp60=Skpeak;
                      Cellrecord{fn}.LaminaPeakTp1Hp60=xx;
                    Cellrecord{fn}.ColumaPeakTp1Hp60=yy;
                    %charge----------------------------------------------------------------------------------
                    avg = gridavg2(X,Y,Zchg,xtick,ytick,false,1,0);
                    xx=nansum(avg);
                    yy=nansum(avg,2);
                    [mxx ii]=max(xx);
                    jj=find(xx>mxx*Perc);
                    
                    if ~isempty(jj)
                        %                      ctext=abs(xi(jj(1))-xi(jj(end)));
                        ctextratio=abs(xi(jj(1))-xi(jj(end)))/abs(xlim(2)-xlim(1));
                        Dx=max(dx,abs(xi(jj(1))-xi(jj(end))));
                    else Dx=0;ctext=0;ctextratio=0;
                    end
                    
                    [myy ii]=max(yy);
                    jj=find(yy>myy*Perc);
                    if ~isempty(jj)
                        Dy=max(dx,abs(yi(jj(1))-yi(jj(end))));
                    else Dy=0;
                    end
             
                    D=[];Skchg=[];Laydist=[];
                    for bi=1:length(Bnd)-1
                        
                             if bi==length(Bnd)-1
                             indBnd=find(xi>=Bnd(bi));
                        else
                            indBnd=find(xi>=Bnd(bi)&xi<Bnd(bi+1));
                        end
                        
                        
                        y1=nansum(avg(:,indBnd),2);
                        s1=skewness(y1);
                        Laydist=[Laydist y1];
                   
                        [mx ii]=max(y1);
                        jj=find(y1>mx*Perc);
                        
                        if ~isempty(jj)
                            D1=max(dx,abs(yi(jj(1))-yi(jj(end))));
                        else D1=0;
                        end
                        D=[D D1];
                        Skchg=[Skchg s1];
                    end
                    
                    DW=[Dx ctextratio  Dy Dy./abs(1000) D];
                    DistWidthCharge=[DistWidthCharge;DW];
                  
                    Skchglayer=[Skchglayer;Skchg];
                    Cellrecord{fn}.LaydistChgTp1Hp60=Laydist;
                    
                    Cellrecord{fn}.DistWidthChargeTp1Hp60=DW;
                    Cellrecord{fn}.SkchglayerTp1Hp60=Skchg;
                    
                    Cellrecord{fn}.LaminaChgTp1Hp60=xx;
                    Cellrecord{fn}.ColumaChgTp1Hp60=yy;
                    %density---------------------------------------------------------------
                    avg = gridavg2(X,Y,Zflag,xtick,ytick,false,1,0);
                    xx=nansum(avg);
                    yy=nansum(avg,2);
                    [mxx ii]=max(xx);
                    jj=find(xx>mxx*Perc);
                    
                    if ~isempty(jj)
                        %                      ctext=abs(xi(jj(1))-xi(jj(end)));
                        ctextratio=abs(xi(jj(1))-xi(jj(end)))/abs(xlim(2)-xlim(1));
                        Dx=max(dx,abs(xi(jj(1))-xi(jj(end))));
                    else Dx=0;ctext=0;ctextratio=0;
                    end
                    [myy ii]=max(yy);
                    jj=find(yy>myy*Perc);
                    if ~isempty(jj)
                        Dy=max(dx,abs(yi(jj(1))-yi(jj(end))));
                    else Dy=0;
                    end
                    
                    D=[];Skdens=[];Laydist=[];
                    for bi=1:length(Bnd)-1
                       
                             if bi==length(Bnd)-1
                             indBnd=find(xi>=Bnd(bi));
                        else
                            indBnd=find(xi>=Bnd(bi)&xi<Bnd(bi+1));
                        end
                        
                        
                        y1=nansum(avg(:,indBnd),2);
                        s1=skewness(y1);
                        Laydist=[Laydist y1];
                        
                        [mx ii]=max(y1);
                        jj=find(y1>mx*Perc);
                        
                        if ~isempty(jj)
                            D1=max(dx,abs(yi(jj(1))-yi(jj(end))));
                        else D1=0;
                        end
                        D=[D D1];
                        Skdens=[Skdens s1];
                    end
                    DW=[Dx ctextratio  Dy Dy./abs(1000) D];
                    DistWidthDense=[DistWidthDense;DW];
                    Cellrecord{fn}.LaydistDensTp1Hp60=Laydist;                 
                    Skdenslayer=[Skdenslayer;Skdens];
                    Cellrecord{fn}.DistWidthdensTp1Hp60=DW;
                    Cellrecord{fn}.SkdenslayerTp1Hp60=Skdens;
                    Cellrecord{fn}.LaminaAreaTp1Hp60=xx;
                    Cellrecord{fn}.ColumaAreaTp1Hp60=yy;
                    
                end
            end
            
            if handles.avgTp==1&&abs(handles.avgHoldingPotential-10)<1&&isfield(Cellrecord{fn},'meanlatencyHighMg10')&&~isempty(Cellrecord{fn}.meanlatencyHighMg10)
                
                Soma= Cellrecord{fn}.SomaCoordinates;
                flipimg=Cellrecord{fn}.flipimg;
                pth=char(Cellrecord{fn}.Pth);
                stimcoordinates=Cellrecord{fn}.StimCoordinates;
                rotate_angle=Cellrecord{fn}.SpatialRotation;
                [stimCoordinates]=Matrixrotate_singlemap_GUI(Soma,stimcoordinates,flipimg,rotate_angle);
                X=stimCoordinates(1,:);
                Y=stimCoordinates(2,:);
                Zpk=zeros(size(X)).*NaN;
                Zchg=zeros(size(X)).*NaN;
                Zlat=zeros(size(X)).*NaN;
                  Ymat=zeros(size(X)).*NaN;
                   Xmat=zeros(size(X)).*NaN;
                if abs(Bnd(indb)-Bnd(1))>=handles.mindist&&abs(Bnd(indb)-Bnd(1))<=handles.maxdist
                    N=N+1;
                    xaxis=[xaxis,X];
                    yaxis=[yaxis,Y];
                    
                    dX=[dX,Cellrecord{fn}.PatternSpacing(1)];
                    Zflag=Cellrecord{fn}.meanflagHighMg10;
                    IndZ=find(Zflag>0);
                    if length(IndZ)
                     dist=sqrt(X(IndZ).^2+Y(IndZ).^2);
                    dist=sort(dist);
                    Dist8010=dist(ceil(length(dist)*0.8));
                    else Dist8010=0;
                    end
                    Zpk(IndZ)=Cellrecord{fn}.meanpeakHighMg10(IndZ);
                    
                    Zchg(IndZ)=Cellrecord{fn}.meanareaHighMg10(IndZ);
                    
                    Zlat(IndZ)=Cellrecord{fn}.meanlatencyHighMg10(IndZ);
                    Ymat(IndZ)=Y(IndZ);
                    Xmat(IndZ)=X(IndZ);
                    Pchg=[];
                    Ppeak=[];
                    Parea=[];
                    chg=[];
                    peak=[];
                    A10=[];
                       Mchg=[];
                    Mpeak=[];
                          Medianchg=[];
                    Medianpeak=[];
                     AxisSum10=[];
                     expSum10=[];
                    ChgexpSum10=[];
                    PkexpSum10=[];
                    
                    Dist8010Layer=[];
                    Xmat=Xmat./sqrt(Xmat.^2+Ymat.^2);
                                Ymat=Ymat./sqrt(Xmat.^2+Ymat.^2);
                    for bi=1:length(Bnd)-1
                        
                             if bi==1
                                 indBnd=find(X<Bnd(bi+1));
                             elseif bi==length(Bnd)-1
                                 indBnd=find(X>=Bnd(bi));
                             else
                            indBnd=find(X>=Bnd(bi)&X<Bnd(bi+1));
                        end
                      
                        
                        Pchg=[Pchg nansum(Zchg(indBnd))/nansum(Zchg)];
                        Ppeak=[Ppeak nansum(Zpk(indBnd))/nansum(Zpk)];
                        Parea=[Parea nansum(Zflag(indBnd))/nansum(Zflag)];
                       chg=[chg nansum(Zchg(indBnd))];
                        peak=[peak nansum(Zpk(indBnd))];
                        Mchg=[Mchg nanmean(Zchg(indBnd))];
                        Mpeak=[Mpeak nanmean(Zpk(indBnd))];
                           Medianchg=[Medianchg nanmedian(Zchg(indBnd))];
                        Medianpeak=[Medianpeak nanmedian(Zpk(indBnd))];
                        A10=[A10 nansum(Zpk(indBnd)>0).*Cellrecord{fn}.PatternSpacing(1).^2];
                         AxisSum10=[AxisSum10 nansum(Ymat(indBnd))];
                         expSum10=[expSum10 abs(nansum((Xmat(indBnd)+i.*Ymat(indBnd))./sqrt(Xmat(indBnd).^2+Ymat(indBnd).^2)))./nansum(abs((Xmat(indBnd)+i.*Ymat(indBnd))./sqrt(Xmat(indBnd).^2+Ymat(indBnd).^2)))];
                        ChgexpSum10=[ChgexpSum10 abs(nansum(Zchg(indBnd).*(Xmat(indBnd)+i.*Ymat(indBnd))./sqrt(Xmat(indBnd).^2+Ymat(indBnd).^2)))./nansum(abs(Zchg(indBnd).*(Xmat(indBnd)+i.*Ymat(indBnd))./sqrt(Xmat(indBnd).^2+Ymat(indBnd).^2)))];
                        PkexpSum10=[PkexpSum10 abs(nansum(Zpk(indBnd).*(Xmat(indBnd)+i.*Ymat(indBnd))./sqrt(Xmat(indBnd).^2+Ymat(indBnd).^2)))./nansum(abs(Zpk(indBnd).*(Xmat(indBnd)+i.*Ymat(indBnd))./sqrt(Xmat(indBnd).^2+Ymat(indBnd).^2)))];
                         indB=find(X(IndZ)>=Bnd(bi)&X(IndZ)<Bnd(bi+1));
                         if length(indB)>0
%                         dist=sqrt(X(IndZ(indB)).^2+Y(IndZ(indB)).^2);
                   dist=abs(Y(IndZ(indB)));
                        dist=sort(dist);
                         %Dist8010Layer=[Dist8010Layer abs(dist(ceil(length(dist)*0.1))-dist(ceil(length(dist)*0.9)))];
                          Dist8010Layer=[Dist8010Layer dist(ceil(length(dist)*0.8))];
                         else
                             Dist8010Layer=[Dist8010Layer 0];
                         end
                    end
                    Pchglayer=[Pchglayer;Pchg];
                    Ppeaklayer=[Ppeaklayer;Ppeak];
                    Parealayer=[Parealayer;Parea];
                    Pk=[Pk,Zpk];
                    Chg=[Chg, Zchg];                    
                    Lat=[Lat,Zlat];
                    Den=[Den Zflag'];
                    Cellrecord{fn}.PchgTp1Hp10=Pchg;
                    Cellrecord{fn}.PpeakTp1Hp10=Ppeak;
                    Cellrecord{fn}.PareaTp1Hp10=Parea;
                    Cellrecord{fn}.chgTp1Hp10=chg;
                    Cellrecord{fn}.peakTp1Hp10=peak;
                    Cellrecord{fn}.AreaTp1Hp10=A10;
                    Cellrecord{fn}.MchgTp1Hp10=Mchg;
                    Cellrecord{fn}.MpeakTp1Hp10=Mpeak;
                    Cellrecord{fn}.AxisSumTp1Hp10=AxisSum10;
                     
                    Cellrecord{fn}.expSum10=expSum10;
                    Cellrecord{fn}.ChgexpSum10=ChgexpSum10;
                    Cellrecord{fn}.PkexpSum10=PkexpSum10;
                             Cellrecord{fn}.MedianchgTp1Hp10=Medianchg;
                    Cellrecord{fn}.MedianpeakTp1Hp10=Medianpeak;
                     Cellrecord{fn}.Dist8010Tp1Hp60=Dist8010;
                     Cellrecord{fn}.Dist8010LayerTp1Hp10=Dist8010Layer;
                    
                    % Peak------------------------------------------
                    dx =Cellrecord{fn}.PatternSpacing(1); dy = dx;
                    xlim = [ceil(min(X))-dx-0.1 ceil(max(X))]; ylim = [ceil(min(Y))-dy-0.1 ceil(max(Y))];
                    
                    
                    xtick = xlim(1):dx:xlim(2); ytick = ylim(1):dy:ylim(2);
                    xi = xtick(1:end-1)+dx/2; yi = ytick(1:end-1)+dy/2;
                    avg = gridavg2(X,Y,Zpk,xtick,ytick,false,1,0);
                                       yy=nansum(avg,2);
                                       xx=nansum(avg);
                    [mxx ii]=max(xx);
                    jj=find(xx>mxx*Perc);
                    
                    if ~isempty(jj)
                        %                      ctext=abs(xi(jj(1))-xi(jj(end)));
                        ctextratio=abs(xi(jj(1))-xi(jj(end)))/abs(xlim(2)-xlim(1));
                        Dx=max(dx,abs(xi(jj(1))-xi(jj(end))));
                    else Dx=0;ctext=0;ctextratio=0;
                    end
                    [myy ii]=max(yy);
                    jj=find(yy>myy*Perc);
                    if ~isempty(jj)
                        Dy=max(dx,abs(yi(jj(1))-yi(jj(end))));
                    else Dy=0;
                    end
                  
                    D=[];Skpeak=[];Laydist=[];
                    for bi=1:length(Bnd)-1
                        
                          if bi==length(Bnd)-1
                             indBnd=find(xi>=Bnd(bi));
                        else
                            indBnd=find(xi>=Bnd(bi)&xi<Bnd(bi+1));
                        end
                        
                        
                        y1=nansum(avg(:,indBnd),2);
                        s1=skewness(y1);
                        Laydist=[Laydist y1];
                     
                        [mx ii]=max(y1);
                        jj=find(y1>mx*Perc);
                        
                        if ~isempty(jj)
                            D1=max(dx,abs(yi(jj(1))-yi(jj(end))));
                        else D1=0;
                        end
                        D=[D D1];
                        Skpeak=[Skpeak s1];
                    end
                   
                    DW=[Dx ctextratio  Dy Dy./abs(ylim(1)-ylim(2)) D];
                    DistWidthPeak=[DistWidthPeak;DW];
                    SkPeaklayer=[SkPeaklayer;Skpeak];
                    Cellrecord{fn}.LaydistPkTp1Hp10=Laydist;
                    Cellrecord{fn}.DistWidthPeakTp1Hp10=DW;
                    Cellrecord{fn}.SkPeaklayerTp1Hp10=Skpeak;
                    
                    Cellrecord{fn}.LaminaPeakTp1Hp10=xx;
                    Cellrecord{fn}.ColumaPeakTp1Hp10=yy;
                    
                     %charge----------------------------------------------------------------------------------
                    avg = gridavg2(X,Y,Zchg,xtick,ytick,false,1,0);
                    xx=nansum(avg);
                    yy=nansum(avg,2);
                    [mxx ii]=max(xx);
                    jj=find(xx>mxx*Perc);
                    
                    if ~isempty(jj)
                        %                      ctext=abs(xi(jj(1))-xi(jj(end)));
                        ctextratio=abs(xi(jj(1))-xi(jj(end)))/abs(xlim(2)-xlim(1));
                        Dx=max(dx,abs(xi(jj(1))-xi(jj(end))));
                    else Dx=0;ctext=0;ctextratio=0;
                    end
                    
                    [myy ii]=max(yy);
                    jj=find(yy>myy*Perc);
                    if ~isempty(jj)
                        Dy=max(dx,abs(yi(jj(1))-yi(jj(end))));
                    else Dy=0;
                    end
             
                    D=[];Skchg=[];Laydist=[];
                    for bi=1:length(Bnd)-1
                       
                         if bi==length(Bnd)-1
                             indBnd=find(xi>=Bnd(bi));
                        else
                            indBnd=find(xi>=Bnd(bi)&xi<Bnd(bi+1));
                        end
                        
                        
                        
                        y1=nansum(avg(:,indBnd),2);
                        s1=skewness(y1);
                        Laydist=[Laydist y1];
                   
                        [mx ii]=max(y1);
                        jj=find(y1>mx*Perc);
                        
                        if ~isempty(jj)
                            D1=max(dx,abs(yi(jj(1))-yi(jj(end))));
                        else D1=0;
                        end
                        D=[D D1];
                        Skchg=[Skchg s1];
                    end
                    
                    DW=[Dx ctextratio  Dy Dy./abs(ylim(1)-ylim(2)) D];
                    DistWidthCharge=[DistWidthCharge;DW];
                  
                    Skchglayer=[Skchglayer;Skchg];
                    Cellrecord{fn}.LaydistChgTp1Hp10=Laydist;                    
                    Cellrecord{fn}.DistWidthChargeTp1Hp10=DW;
                    Cellrecord{fn}.SkchglayerTp1Hp10=Skchg;
                      
                    Cellrecord{fn}.LaminaChgTp1Hp10=xx;
                    Cellrecord{fn}.ColumaChgTp1Hp10=yy;
                    
                    %density---------------------------------------------------------------
                    avg = gridavg2(X,Y,Zflag,xtick,ytick,false,1,0);
                    xx=nansum(avg);
                    yy=nansum(avg,2);
                    [mxx ii]=max(xx);
                    jj=find(xx>mxx*Perc);
                    
                    if ~isempty(jj)
                        %                      ctext=abs(xi(jj(1))-xi(jj(end)));
                        ctextratio=abs(xi(jj(1))-xi(jj(end)))/abs(xlim(2)-xlim(1));
                        Dx=max(dx,abs(xi(jj(1))-xi(jj(end))));
                    else Dx=0;ctext=0;ctextratio=0;
                    end
                    [myy ii]=max(yy);
                    jj=find(yy>myy*Perc);
                    if ~isempty(jj)
                        Dy=max(dx,abs(yi(jj(1))-yi(jj(end))));
                    else Dy=0;
                    end
                    
                    D=[];Skdens=[];Laydist=[];
                    for bi=1:length(Bnd)-1
                        
                            if bi==length(Bnd)-1
                             indBnd=find(xi>=Bnd(bi));
                        else
                            indBnd=find(xi>=Bnd(bi)&xi<Bnd(bi+1));
                        end
                        
                     
                        
                        y1=nansum(avg(:,indBnd),2);
                        s1=skewness(y1);
                        Laydist=[Laydist y1];
                        
                        [mx ii]=max(y1);
                        jj=find(y1>mx*Perc);
                        
                        if ~isempty(jj)
                            D1=max(dx,abs(yi(jj(1))-yi(jj(end))));
                        else D1=0;
                        end
                        D=[D D1];
                        Skdens=[Skdens s1];
                    end
                    DW=[Dx ctextratio  Dy Dy./abs(ylim(1)-ylim(2)) D];
                    DistWidthDense=[DistWidthDense;DW];
                    Cellrecord{fn}.LaydistDensTp1Hp10=Laydist;                 
                    Skdenslayer=[Skdenslayer;Skdens];
                    Cellrecord{fn}.DistWidthdensTp1Hp10=DW;
                    Cellrecord{fn}.SkdenslayerTp1Hp10=Skdens;
                      
                    Cellrecord{fn}.LaminaAreaTp1Hp10=xx;
                    Cellrecord{fn}.ColumaAreaTp1Hp10=yy;
                end
            end
            
            
            
            if handles.avgTp==2&&abs(handles.avgHoldingPotential+60)<1&&isfield(Cellrecord{fn},'meanlatencyPtxInd60')&&~isempty(Cellrecord{fn}.meanlatencyPtxInd60)
                
                Soma= Cellrecord{fn}.SomaCoordinates;
                flipimg=Cellrecord{fn}.flipimg;
                pth=char(Cellrecord{fn}.Pth);
                stimcoordinates=Cellrecord{fn}.StimCoordinates;
                rotate_angle=Cellrecord{fn}.SpatialRotation;
                [stimCoordinates]=Matrixrotate_singlemap_GUI(Soma,stimcoordinates,flipimg,rotate_angle);
                X=stimCoordinates(1,:);
                Y=stimCoordinates(2,:);
                Zpk=zeros(size(X)).*NaN;
                Zchg=zeros(size(X)).*NaN;
                Zlat=zeros(size(X)).*NaN;
                
                if abs(Bnd(indb)-Bnd(1))>=handles.mindist&&abs(Bnd(indb)-Bnd(1))<=handles.maxdist
                    
                    N=N+1;
                    xaxis=[xaxis,X];
                    yaxis=[yaxis,Y];
                    dX=[dX,Cellrecord{fn}.PatternSpacing(1)];
                    Zflag=Cellrecord{fn}.meanflagPtxInd60;
                    IndZ=find(Zflag>0);
                    
                    Zpk(IndZ)=Cellrecord{fn}.meanpeakPtxInd60(IndZ);
                    Zchg(IndZ)=Cellrecord{fn}.meanareaPtxInd60(IndZ);
                    Zlat(IndZ)=Cellrecord{fn}.meanlatencyPtxInd60(IndZ);
                    
                    
                    Pchg=[];
                    Ppeak=[];
                    Parea=[];
                    chg=[];
                    peak=[];
                    A60=[];
                     Mchg=[];
                    Mpeak=[];
                    for bi=1:length(Bnd)-1
                        
                            indBnd=find(X>=Bnd(bi)&X<Bnd(bi+1));
                      
                        
                        Pchg=[Pchg nansum(Zchg(indBnd))/nansum(Zchg)];
                        Ppeak=[Ppeak nansum(Zpk(indBnd))/nansum(Zpk)];
                        Parea=[Parea nansum(Zflag(indBnd))/nansum(Zflag)];
                        chg=[chg nansum(Zchg(indBnd))];
                        peak=[peak nansum(Zpk(indBnd))];
                        Mchg=[Mchg nanmean(Zchg(indBnd))];
                        Mpeak=[Mpeak nanmean(Zpk(indBnd))];
                        A60=[A60 nansum(Zpk(indBnd)>0).*Cellrecord{fn}.PatternSpacing(1).^2];
                    end
                    Pk=[Pk,Zpk];
                    Pchglayer=[Pchglayer;Pchg];
                    Ppeaklayer=[Ppeaklayer;Ppeak];
                     Parealayer=[Parealayer;Parea];
                    Chg=[Chg, Zchg];
                    Lat=[Lat,Zlat];
                    Den=[Den Zflag'];
                    Cellrecord{fn}.PchgTp2Hp60=Pchg;
                    Cellrecord{fn}.PpeakTp2Hp60=Ppeak;
                    Cellrecord{fn}.PareaTp2Hp60=Parea;
                    Cellrecord{fn}.chgTp2Hp60=chg;
                    Cellrecord{fn}.peakTp2Hp60=peak;
                    Cellrecord{fn}.AreaTp2Hp60=A60;
                     Cellrecord{fn}.MchgTp2Hp60=Mchg;
                    Cellrecord{fn}.MpeakTp2Hp60=Mpeak;
                    
                    ZflagD=Cellrecord{fn}.meanflagDPtxInd60;
                    IndZD=find(ZflagD>0); 
                    Cellrecord{fn}.DirectAreaTp2Hp60=length(IndZD).*Cellrecord{fn}.PatternSpacing(1).^2;
                    DirArea=[DirArea;length(IndZD).*Cellrecord{fn}.PatternSpacing(1).^2];
                    % Peak------------------------------------------
                    dx =Cellrecord{fn}.PatternSpacing(1); dy = dx;
                    xlim = [ceil(min(X))-dx-0.1 ceil(max(X))]; ylim = [ceil(min(Y))-dy-0.1 ceil(max(Y))];
                    
                    
                    xtick = xlim(1):dx:xlim(2); ytick = ylim(1):dy:ylim(2);
                    xi = xtick(1:end-1)+dx/2; yi = ytick(1:end-1)+dy/2;
                    avg = gridavg2(X,Y,Zpk,xtick,ytick,false,1,0);
                                       yy=nansum(avg,2);
                                       xx=nansum(avg);
                    [mxx ii]=max(xx);
                    jj=find(xx>mxx*Perc);
                    
                    if ~isempty(jj)
                        %                      ctext=abs(xi(jj(1))-xi(jj(end)));
                        ctextratio=abs(xi(jj(1))-xi(jj(end)))/abs(xlim(2)-xlim(1));
                        Dx=max(dx,abs(xi(jj(1))-xi(jj(end))));
                    else Dx=0;ctext=0;ctextratio=0;
                    end
                    [myy ii]=max(yy);
                    jj=find(yy>myy*Perc);
                    if ~isempty(jj)
                        Dy=max(dx,abs(yi(jj(1))-yi(jj(end))));
                    else Dy=0;
                    end
                  
                    D=[];Skpeak=[];Laydist=[];
                    for bi=1:length(Bnd)-1
                        
                          if bi==length(Bnd)-1
                             indBnd=find(xi>=Bnd(bi));
                        else
                            indBnd=find(xi>=Bnd(bi)&xi<Bnd(bi+1));
                        end
                        
                       
                        
                        y1=nansum(avg(:,indBnd),2);
                        s1=skewness(y1);
                        Laydist=[Laydist y1];
                     
                        [mx ii]=max(y1);
                        jj=find(y1>mx*Perc);
                        
                        if ~isempty(jj)
                            D1=max(dx,abs(yi(jj(1))-yi(jj(end))));
                        else D1=0;
                        end
                        D=[D D1];
                        Skpeak=[Skpeak s1];
                    end
                   
                    DW=[Dx ctextratio  Dy Dy./abs(ylim(1)-ylim(2)) D];
                    DistWidthPeak=[DistWidthPeak;DW];
                    SkPeaklayer=[SkPeaklayer;Skpeak];
                    Cellrecord{fn}.LaydistPkTp2Hp60=Laydist;
                    Cellrecord{fn}.DistWidthPeakTp2Hp60=DW;
                    Cellrecord{fn}.SkPeaklayerTp2Hp60=Skpeak;
                    
                    
                    
                     %charge----------------------------------------------------------------------------------
                    avg = gridavg2(X,Y,Zchg,xtick,ytick,false,1,0);
                    xx=nansum(avg);
                    yy=nansum(avg,2);
                    [mxx ii]=max(xx);
                    jj=find(xx>mxx*Perc);
                    
                    if ~isempty(jj)
                        %                      ctext=abs(xi(jj(1))-xi(jj(end)));
                        ctextratio=abs(xi(jj(1))-xi(jj(end)))/abs(xlim(2)-xlim(1));
                        Dx=max(dx,abs(xi(jj(1))-xi(jj(end))));
                    else Dx=0;ctext=0;ctextratio=0;
                    end
                    
                    [myy ii]=max(yy);
                    jj=find(yy>myy*Perc);
                    if ~isempty(jj)
                        Dy=max(dx,abs(yi(jj(1))-yi(jj(end))));
                    else Dy=0;
                    end
             
                    D=[];Skchg=[];Laydist=[];
                    for bi=1:length(Bnd)-1
                        
                           if bi==length(Bnd)-1
                             indBnd=find(xi>=Bnd(bi));
                        else
                            indBnd=find(xi>=Bnd(bi)&xi<Bnd(bi+1));
                        end
                        
                     
                        
                        y1=nansum(avg(:,indBnd),2);
                        s1=skewness(y1);
                        Laydist=[Laydist y1];
                   
                        [mx ii]=max(y1);
                        jj=find(y1>mx*Perc);
                        
                        if ~isempty(jj)
                            D1=max(dx,abs(yi(jj(1))-yi(jj(end))));
                        else D1=0;
                        end
                        D=[D D1];
                        Skchg=[Skchg s1];
                    end
                    
                    DW=[Dx ctextratio  Dy Dy./abs(ylim(1)-ylim(2)) D];
                    DistWidthCharge=[DistWidthCharge;DW];
                  
                    Skchglayer=[Skchglayer;Skchg];
                    Cellrecord{fn}.LaydistChgTp2Hp60=Laydist;                    
                    Cellrecord{fn}.DistWidthChargeTp2Hp60=DW;
                    Cellrecord{fn}.SkchglayerTp2Hp60=Skchg;
                    
                    
                    %density---------------------------------------------------------------
                    avg = gridavg2(X,Y,Zflag,xtick,ytick,false,1,0);
                    xx=nansum(avg);
                    yy=nansum(avg,2);
                    [mxx ii]=max(xx);
                    jj=find(xx>mxx*Perc);
                    
                    if ~isempty(jj)
                        %                      ctext=abs(xi(jj(1))-xi(jj(end)));
                        ctextratio=abs(xi(jj(1))-xi(jj(end)))/abs(xlim(2)-xlim(1));
                        Dx=max(dx,abs(xi(jj(1))-xi(jj(end))));
                    else Dx=0;ctext=0;ctextratio=0;
                    end
                    [myy ii]=max(yy);
                    jj=find(yy>myy*Perc);
                    if ~isempty(jj)
                        Dy=max(dx,abs(yi(jj(1))-yi(jj(end))));
                    else Dy=0;
                    end
                    
                    D=[];Skdens=[];Laydist=[];
                    for bi=1:length(Bnd)-1
                      
                          if bi==length(Bnd)-1
                             indBnd=find(xi>=Bnd(bi));
                        else
                            indBnd=find(xi>=Bnd(bi)&xi<Bnd(bi+1));
                        end
                        
                       
                        
                        y1=nansum(avg(:,indBnd),2);
                        s1=skewness(y1);
                        Laydist=[Laydist y1];
                        
                        [mx ii]=max(y1);
                        jj=find(y1>mx*Perc);
                        
                        if ~isempty(jj)
                            D1=max(dx,abs(yi(jj(1))-yi(jj(end))));
                        else D1=0;
                        end
                        D=[D D1];
                        Skdens=[Skdens s1];
                    end
                    DW=[Dx ctextratio  Dy Dy./abs(ylim(1)-ylim(2)) D];
                    DistWidthDense=[DistWidthDense;DW];
                    Cellrecord{fn}.LaydistDensTp2Hp60=Laydist;                 
                    Skdenslayer=[Skdenslayer;Skdens];
                    Cellrecord{fn}.DistWidthdensTp2Hp60=DW;
                    Cellrecord{fn}.SkdenslayerTp2Hp60=Skdens;
                    
                end
            end
            
            
            
            
             if handles.avgTp==2&&isfield(Cellrecord{fn},'meanpeakPtxInd50')&&isfield(Cellrecord{fn},'meanpeakPtxInd60')
                 Soma= Cellrecord{fn}.SomaCoordinates;
                flipimg=Cellrecord{fn}.flipimg;
                pth=char(Cellrecord{fn}.Pth);
                stimcoordinates=Cellrecord{fn}.StimCoordinates;
                rotate_angle=Cellrecord{fn}.SpatialRotation;
                [stimCoordinates]=Matrixrotate_singlemap_GUI(Soma,stimcoordinates,flipimg,rotate_angle);
                X=stimCoordinates(1,:);
                Y=stimCoordinates(2,:);
                Zpk=zeros(size(X)).*NaN;
                Zchg=zeros(size(X)).*NaN;
                Zlat=zeros(size(X)).*NaN;
                Zflag=Cellrecord{fn}.Matrix_NMDAonly;
                IndZ=find(Zflag>0);
                if abs(Bnd(indb)-Bnd(1))>=handles.mindist&&abs(Bnd(indb)-Bnd(1))<=handles.maxdist
                  
                      
                    Nsi=Nsi+1;
                 xaxissi=[xaxissi,X];
                yaxissi=[yaxissi,Y];
                AGEsi=[AGEsi,Cellrecord{fn}.age];
                BDsi=[BDsi;Bnd];
                dXsi=[dXsi,Cellrecord{fn}.PatternSpacing(1)];
               
                
               
                Zpk(IndZ)=Cellrecord{fn}.meanpeakPtxInd50(IndZ);
                
                Zchg(IndZ)=Cellrecord{fn}.meanareaPtxInd50(IndZ);
                
                Zlat(IndZ)=Cellrecord{fn}.meanlatencyPtxInd50(IndZ);
                
                ind1=find(X>=Bnd(1)&X<Bnd(2));
                 ind23=find(X>=Bnd(2)&X<Bnd(3));
                  ind4=find(X>=Bnd(3)&X<Bnd(4));
                 ind56=find(X>=Bnd(4)&X<Bnd(5));
                 if ~isempty(IndZ)
                Pchg=[nansum(Zchg(ind1))/nansum(Zchg) nansum(Zchg(ind23))/nansum(Zchg) nansum(Zchg(ind4))/nansum(Zchg) nansum(Zchg(ind56))/nansum(Zchg)];
                Ppeak=[nansum(Zpk(ind1))/nansum(Zpk) nansum(Zpk(ind23))/nansum(Zpk) nansum(Zpk(ind4))/nansum(Zpk) nansum(Zpk(ind56))/nansum(Zpk)];
                 Parea=[nansum(Zflag(ind1))/nansum(Zflag) nansum(Zflag(ind23))/nansum(Zflag) nansum(Zflag(ind4))/nansum(Zflag) nansum(Zflag(ind56))/nansum(Zflag)];
                  chg=[nansum(Zchg(ind1)) nansum(Zchg(ind23)) nansum(Zchg(ind4)) nansum(Zchg(ind56))];
                peak=[nansum(Zpk(ind1)) nansum(Zpk(ind23)) nansum(Zpk(ind4)) nansum(Zpk(ind56))];
                Mchg=[nanmean(Zchg(ind1)) nanmean(Zchg(ind23)) nanmean(Zchg(ind4)) nanmean(Zchg(ind56))];
                Mpeak=[nanmean(Zpk(ind1)) nanmean(Zpk(ind23)) nanmean(Zpk(ind4)) nanmean(Zpk(ind56))];
                Asi=[nansum(Zpk(ind1)>0) nansum(Zpk(ind23)>0) nansum(Zpk(ind4)>0) nansum(Zpk(ind56)>0)].*Cellrecord{fn}.PatternSpacing(1).^2;
                 else
                     Pchg=[NaN NaN NaN NaN];
                     Ppeak=[NaN NaN NaN NaN];
                      Parea=[NaN NaN NaN NaN];
                     chg=[NaN NaN NaN NaN];
                     peak=[NaN NaN NaN NaN];
                     Asi=[NaN NaN NaN NaN];
                      Mchg=[NaN NaN NaN NaN];
                     Mpeak=[NaN NaN NaN NaN];
                     
                 end
                Pksi=[Pksi,Zpk];
                  Pchgsilayer=[Pchgsilayer;Pchg];
                Ppeaksilayer=[Ppeaksilayer;Ppeak];
                 Pareasilayer=[Pareasilayer;Parea];
                
                Chgsi=[Chgsi, Zchg];
                
                Latsi=[Latsi,Zlat];
                Densi=[Densi Zflag'];
                Cellrecord{fn}.PchgTp2si=Pchg;
                      Cellrecord{fn}.PpeakTp2si=Ppeak;
                      Cellrecord{fn}.PareaTp2si=Parea;
                      Cellrecord{fn}.chgTp2si=chg;
                      Cellrecord{fn}.peakTp2si=peak;
                      Cellrecord{fn}.AreaTp2si=Asi;
                  Cellrecord{fn}.MchgTp2si=Mchg;
                      Cellrecord{fn}.MpeakTp2si=Mpeak;
                      
                      dx =Cellrecord{fn}.PatternSpacing(1); dy = dx;
                    xlim = [ceil(min(X))-dx-0.1 ceil(max(X))]; ylim = [ceil(min(Y))-dy-0.1 ceil(max(Y))];
                    
                    
                    xtick = xlim(1):dx:xlim(2); ytick = ylim(1):dy:ylim(2);
                    xi = xtick(1:end-1)+dx/2; yi = ytick(1:end-1)+dy/2;
                    %pk
                     avg = gridavg2(X,Y,Zpk,xtick,ytick,false,1,0); 
                    xx=nansum(avg);
                    yy=nansum(avg,2);
                    [mxx ii]=max(xx);
                    jj=find(xx>mxx*Perc);
                  if ~isempty(jj)
                     ctext=abs(xi(jj(1))-Bnd(4));
                     ctextratio=abs(xi(jj(1))-Bnd(4))/abs(Bnd(1)-Bnd(4));
                      Dx=max(dx,abs(xi(jj(1))-xi(jj(end))));
                 else Dx=0;ctext=0;ctextratio=0;
                 end
                [myy ii]=max(yy);
                    jj=find(yy>myy*Perc);
                 if ~isempty(jj)
                      Dy=max(dx,abs(yi(jj(1))-yi(jj(end))));
                 else Dy=0;
                 end
                     ind1=find(xi>=Bnd(1)&xi<Bnd(2));
                     ind23=find(xi>=Bnd(2)&xi<Bnd(3));
                  ind4=find(xi>=Bnd(3)&xi<Bnd(4));
                 ind56=find(xi>=Bnd(4)&xi<Bnd(5));
                 y1=nansum(avg(:,ind1),2);
%                  [mx ii]=max(y1);
                
%                  jj=find(y1(1:ii-1)<=mx*Perc);
%                  if isempty(jj)
%                      jend=1;
%                  else jend=jj(end);
%                  end
%                  jj=find(y1(ii+1:end)<=mx*Perc);
%                  if isempty(jj)
%                      j1=length(y1);
%                  else j1=jj(1);
%                  end
%                  D1=abs(yi(j1)-yi(jend));
                [mx ii]=max(y1);
                jj=find(y1>mx*Perc);
                 if ~isempty(jj)
                      D1=max(dx,abs(yi(jj(1))-yi(jj(end))));
                 else D1=0;
                 end
                  y23=nansum(avg(:,ind23),2)/length(ind23);
%                  [mx ii]=max(y23);
%                  jj=find(y23(1:ii-1)<=mx*Perc);
%                  if isempty(jj)
%                      jend=1;
%                  else jend=jj(end);
%                  end
%                  jj=find(y23(ii+1:end)<=mx*Perc);
%                  if isempty(jj)
%                      j1=length(y23);
%                  else j1=jj(1);
%                  end
%                  D23=abs(yi(j1)-yi(jend));
                 [mx ii]=max(y23);
                jj=find(y23>mx*Perc);
                 if ~isempty(jj)
                      D23=max(dx,abs(yi(jj(1))-yi(jj(end))));
                 else D23=0;
                 end
                    y4=nansum(avg(:,ind4),2)/length(ind4);
%                  [mx ii]=max(y4);
%                  jj=find(y4(1:ii-1)<=mx*Perc);
%                  if isempty(jj)
%                      jend=1;
%                  else jend=jj(end);
%                  end
%                  jj=find(y4(ii+1:end)<=mx*Perc);
%                  if isempty(jj)
%                      j1=length(y4);
%                  else j1=jj(1);
%                  end
%                  D4=abs(yi(j1)-yi(jend));
                  [mx ii]=max(y4);
                jj=find(y4>mx*Perc);
                 if ~isempty(jj)
                      D4=max(dx,abs(yi(jj(1))-yi(jj(end))));
                 else D4=0;
                 end
                    y56=nansum(avg(:,ind56),2)/length(ind56);
%                  [mx ii]=max(y56);
%                  jj=find(y56(1:ii-1)<=mx*Perc);
%                  if isempty(jj)
%                      jend=1;
%                  else jend=jj(end);
%                  end
%                  jj=find(y56(ii+1:end)<=mx*Perc);
%                  if isempty(jj)
%                      j1=length(y56);
%                  else j1=jj(1);
%                  end
%                  D56=abs(yi(j1)-yi(jend));
                 [mx ii]=max(y56);
                jj=find(y56>mx*Perc);
                 if ~isempty(jj)
                      D56=max(dx,abs(yi(jj(1))-yi(jj(end))));
                 else D56=0;
                 end
              DW=[Dx Dx./abs(Bnd(end)-Bnd(1)) ctext ctextratio  Dy Dy./abs(ylim(1)-ylim(2)) D1 D23 D4 D56];
                 DistWidthPeaksi=[DistWidthPeaksi;DW];
%                  s1=skewness(y1);
%                   s23=skewness(y23);
%                    s4=skewness(y4);
%                     s56=skewness(y56);
%                     Skpeak=[s1 s23 s4 s56];
%                  SkPeaklayersi=[SkPeaklayersi;Skpeak];
                 Cellrecord{fn}.L1PkTp2Hp50si=y1;
                 Cellrecord{fn}.L23PkTp2Hp50si=y23;
                 Cellrecord{fn}.L4PkTp2Hp50si=y4;
                 Cellrecord{fn}.L56PkTp2Hp50si=y56;
                   Cellrecord{fn}.DistWidthPeakTp2Hp50si=DW;
                  
%                    Cellrecord{fn}.SkPeaklayerTp2Hp50si=Skpeak;
                 
                 %charge
                 avg = gridavg2(X,Y,Zchg,xtick,ytick,false,1,0); 
                          xx=nansum(avg);
                    yy=nansum(avg,2);
                    [mxx ii]=max(xx);
                    jj=find(xx>mxx*Perc);
                if ~isempty(jj)
                     ctext=abs(xi(jj(1))-Bnd(4));
                     ctextratio=abs(xi(jj(1))-Bnd(4))/abs(Bnd(1)-Bnd(4));
                      Dx=max(dx,abs(xi(jj(1))-xi(jj(end))));
                 else Dx=0;ctext=0;ctextratio=0;
                 end
                [myy ii]=max(yy);
                    jj=find(yy>myy*Perc);
                 if ~isempty(jj)
                      Dy=max(dx,abs(yi(jj(1))-yi(jj(end))));
                 else Dy=0;
                 end
                     ind1=find(xi>=Bnd(1)&xi<Bnd(2));
                     ind23=find(xi>=Bnd(2)&xi<Bnd(3));
                  ind4=find(xi>=Bnd(3)&xi<Bnd(4));
                 ind56=find(xi>=Bnd(4)&xi<Bnd(5));
                 y1=nansum(avg(:,ind1),2);
%                  [mx ii]=max(y1);
                
%                  jj=find(y1(1:ii-1)<=mx*Perc);
%                  if isempty(jj)
%                      jend=1;
%                  else jend=jj(end);
%                  end
%                  jj=find(y1(ii+1:end)<=mx*Perc);
%                  if isempty(jj)
%                      j1=length(y1);
%                  else j1=jj(1);
%                  end
%                  D1=abs(yi(j1)-yi(jend));
                [mx ii]=max(y1);
                jj=find(y1>mx*Perc);
                 if ~isempty(jj)
                      D1=max(dx,abs(yi(jj(1))-yi(jj(end))));
                 else D1=0;
                 end
                  y23=nansum(avg(:,ind23),2);
%                  [mx ii]=max(y23);
%                  jj=find(y23(1:ii-1)<=mx*Perc);
%                  if isempty(jj)
%                      jend=1;
%                  else jend=jj(end);
%                  end
%                  jj=find(y23(ii+1:end)<=mx*Perc);
%                  if isempty(jj)
%                      j1=length(y23);
%                  else j1=jj(1);
%                  end
%                  D23=abs(yi(j1)-yi(jend));
                 [mx ii]=max(y23);
                jj=find(y23>mx*Perc);
                 if ~isempty(jj)
                      D23=max(dx,abs(yi(jj(1))-yi(jj(end))));
                 else D23=0;
                 end
                 %length(ind4);
%                  [mx ii]=max(y4);
%                  jj=find(y4(1:ii-1)<=mx*Perc);
%                  if isempty(jj)
%                      jend=1;
%                  else jend=jj(end);
%                  end
%                  jj=find(y4(ii+1:end)<=mx*Perc);
%                  if isempty(jj)
%                      j1=length(y4);
%                  else j1=jj(1);
%                  end
%                  D4=abs(yi(j1)-yi(jend));
                  [mx ii]=max(y4);
                jj=find(y4>mx*Perc);
                 if ~isempty(jj)
                      D4=max(dx,abs(yi(jj(1))-yi(jj(end))));
                 else D4=0;
                 end
                    y56=nansum(avg(:,ind56),2);
%                  [mx ii]=max(y56);
%                  jj=find(y56(1:ii-1)<=mx*Perc);
%                  if isempty(jj)
%                      jend=1;
%                  else jend=jj(end);
%                  end
%                  jj=find(y56(ii+1:end)<=mx*Perc);
%                  if isempty(jj)
%                      j1=length(y56);
%                  else j1=jj(1);
%                  end
%                  D56=abs(yi(j1)-yi(jend));
                 [mx ii]=max(y56);
                jj=find(y56>mx*Perc);
                 if ~isempty(jj)
                      D56=max(dx,abs(yi(jj(1))-yi(jj(end))));
                 else D56=0;
                 end
              DW=[Dx Dx./abs(Bnd(end)-Bnd(1)) ctext ctextratio  Dy Dy./abs(ylim(1)-ylim(2)) D1 D23 D4 D56];
                 DistWidthChargesi=[DistWidthChargesi;DW];
%                  s1=skewness(y1);
%                   s23=skewness(y23);
%                    s4=skewness(y4);
%                     s56=skewness(y56);
%                     Skchg=[s1 s23 s4 s56];
%                  Skchglayer=[Skchglayer;Skchg];
                 Cellrecord{fn}.L1chgTp2Hp50si=y1;
                 Cellrecord{fn}.L23chgTp2Hp50si=y23;
                 Cellrecord{fn}.L4chgTp2Hp50si=y4;
                 Cellrecord{fn}.L56chgTp2Hp50si=y56;
                 Cellrecord{fn}.DistWidthChargeTp2Hp50si=DW;
%                    Cellrecord{fn}.SkchglayerTp2Hp50si=Skchg;
                 %density
                  avg = gridavg2(X,Y,Zflag,xtick,ytick,false,1,0);
                  
                                   xx=nansum(avg);
                    yy=nansum(avg,2);
                    [mxx ii]=max(xx);
                    jj=find(xx>mxx*Perc);
                  if ~isempty(jj)
                     ctext=abs(xi(jj(1))-Bnd(4));
                     ctextratio=abs(xi(jj(1))-Bnd(4))/abs(Bnd(1)-Bnd(4));
                      Dx=max(dx,abs(xi(jj(1))-xi(jj(end))));
                 else Dx=0;ctext=0;ctextratio=0;
                 end
                [myy ii]=max(yy);
                    jj=find(yy>myy*Perc);
                 if ~isempty(jj)
                      Dy=max(dx,abs(yi(jj(1))-yi(jj(end))));
                 else Dy=0;
                 end
                     ind1=find(xi>=Bnd(1)&xi<Bnd(2));
                     ind23=find(xi>=Bnd(2)&xi<Bnd(3));
                  ind4=find(xi>=Bnd(3)&xi<Bnd(4));
                 ind56=find(xi>=Bnd(4)&xi<Bnd(5));
                 y1=nansum(avg(:,ind1),2);
%                  [mx ii]=max(y1);
                
%                  jj=find(y1(1:ii-1)<=mx*Perc);
%                  if isempty(jj)
%                      jend=1;
%                  else jend=jj(end);
%                  end
%                  jj=find(y1(ii+1:end)<=mx*Perc);
%                  if isempty(jj)
%                      j1=length(y1);
%                  else j1=jj(1);
%                  end
%                  D1=abs(yi(j1)-yi(jend));
                [mx ii]=max(y1);
                jj=find(y1>mx*Perc);
                 if ~isempty(jj)
                      D1=max(dx,abs(yi(jj(1))-yi(jj(end))));
                 else D1=0;
                 end
                  y23=nansum(avg(:,ind23),2);
%                  [mx ii]=max(y23);
%                  jj=find(y23(1:ii-1)<=mx*Perc);
%                  if isempty(jj)
%                      jend=1;
%                  else jend=jj(end);
%                  end
%                  jj=find(y23(ii+1:end)<=mx*Perc);
%                  if isempty(jj)
%                      j1=length(y23);
%                  else j1=jj(1);
%                  end
%                  D23=abs(yi(j1)-yi(jend));
                 [mx ii]=max(y23);
                jj=find(y23>mx*Perc);
                 if ~isempty(jj)
                      D23=max(dx,abs(yi(jj(1))-yi(jj(end))));
                 else D23=0;
                 end
                    y4=nansum(avg(:,ind4),2);
%                  [mx ii]=max(y4);
%                  jj=find(y4(1:ii-1)<=mx*Perc);
%                  if isempty(jj)
%                      jend=1;
%                  else jend=jj(end);
%                  end
%                  jj=find(y4(ii+1:end)<=mx*Perc);
%                  if isempty(jj)
%                      j1=length(y4);
%                  else j1=jj(1);
%                  end
%                  D4=abs(yi(j1)-yi(jend));
                  [mx ii]=max(y4);
                jj=find(y4>mx*Perc);
                 if ~isempty(jj)
                      D4=max(dx,abs(yi(jj(1))-yi(jj(end))));
                 else D4=0;
                 end
                    y56=nansum(avg(:,ind56),2);
%                  [mx ii]=max(y56);
%                  jj=find(y56(1:ii-1)<=mx*Perc);
%                  if isempty(jj)
%                      jend=1;
%                  else jend=jj(end);
%                  end
%                  jj=find(y56(ii+1:end)<=mx*Perc);
%                  if isempty(jj)
%                      j1=length(y56);
%                  else j1=jj(1);
%                  end
%                  D56=abs(yi(j1)-yi(jend));
                 [mx ii]=max(y56);
                jj=find(y56>mx*Perc);
                 if ~isempty(jj)
                      D56=max(dx,abs(yi(jj(1))-yi(jj(end))));
                 else D56=0;
                 end
              DW=[Dx Dx./abs(Bnd(end)-Bnd(1)) ctext ctextratio  Dy Dy./abs(ylim(1)-ylim(2)) D1 D23 D4 D56];
                
                 DistWidthdenssi=[DistWidthdenssi;DW];
%                  s1=skewness(y1);
%                   s23=skewness(y23);
%                    s4=skewness(y4);
%                     s56=skewness(y56);
%                     Skdens=[s1 s23 s4 s56];
%                  Skdenslayer=[Skdenslayer;Skdens];
                  Cellrecord{fn}.L1densTp2Hp50si=y1;
                 Cellrecord{fn}.L23densTp2Hp50si=y23;
                 Cellrecord{fn}.L4densTp2Hp50si=y4;
                 Cellrecord{fn}.L56densTp2Hp50si=y56;
                 Cellrecord{fn}.DistWidthdensTp2Hp50si=DW;
%                    Cellrecord{fn}.SkdenslayerTp2Hp50=Skdens;
                 
                end
            end
                
                
                
                
                
                
                
                
            if handles.avgTp==2&&abs(handles.avgHoldingPotential-50)<1&&isfield(Cellrecord{fn},'meanlatencyPtxInd50')&&~isempty(Cellrecord{fn}.meanlatencyPtxInd50)
            
                Soma= Cellrecord{fn}.SomaCoordinates;
                flipimg=Cellrecord{fn}.flipimg;
                pth=char(Cellrecord{fn}.Pth);
                stimcoordinates=Cellrecord{fn}.StimCoordinates;
                rotate_angle=Cellrecord{fn}.SpatialRotation;
                [stimCoordinates]=Matrixrotate_singlemap_GUI(Soma,stimcoordinates,flipimg,rotate_angle);
                X=stimCoordinates(1,:);
                Y=stimCoordinates(2,:);
                Zpk=zeros(size(X)).*NaN;
                Zchg=zeros(size(X)).*NaN;
                Zlat=zeros(size(X)).*NaN;
               
                   
                if abs(Bnd(indb)-Bnd(1))>=handles.mindist&&abs(Bnd(indb)-Bnd(1))<=handles.maxdist
                     N=N+1;
               
                xaxis=[xaxis,X];
                yaxis=[yaxis,Y];
               
                dX=[dX,Cellrecord{fn}.PatternSpacing(1)];
                Zflag=Cellrecord{fn}.meanflagPtxInd50;
                IndZ=find(Zflag>0);
                Zpk(IndZ)=Cellrecord{fn}.meanpeakPtxInd50(IndZ);
                
                Zchg(IndZ)=Cellrecord{fn}.meanareaPtxInd50(IndZ);
                
                Zlat(IndZ)=Cellrecord{fn}.meanlatencyPtxInd50(IndZ);
                
                  Pchg=[];
                    Ppeak=[];
                    Parea=[];
                    chg=[];
                    peak=[];
                    A50=[];
                     Mchg=[];
                    Mpeak=[];
                    
                    for bi=1:length(Bnd)-1
                        
                            indBnd=find(X>=Bnd(bi)&X<Bnd(bi+1));
                       
                        Pchg=[Pchg nansum(Zchg(indBnd))/nansum(Zchg)];
                        Ppeak=[Ppeak nansum(Zpk(indBnd))/nansum(Zpk)];
                           Parea=[Parea nansum(Zflag(indBnd))/nansum(Zflag)];
                        Mchg=[Mchg nanmean(Zchg(indBnd))];
                        Mpeak=[Mpeak nanmean(Zpk(indBnd))];
                        A50=[A50 nansum(Zpk(indBnd)>0).*Cellrecord{fn}.PatternSpacing(1).^2];
                          chg=[chg nansum(Zchg(indBnd))];
                        peak=[peak nansum(Zpk(indBnd))];
                    end
                    
                Pk=[Pk,Zpk];
                  Pchglayer=[Pchglayer;Pchg];
                Ppeaklayer=[Ppeaklayer;Ppeak];
                 Parealayer=[Parealayer;Parea];
                Chg=[Chg, Zchg];
                
                Lat=[Lat,Zlat];
                Den=[Den Zflag'];
                Cellrecord{fn}.PchgTp2Hp50=Pchg;
                      Cellrecord{fn}.PpeakTp2Hp50=Ppeak;
                      Cellrecord{fn}.AreaTp2Hp50=A50;
                      Cellrecord{fn}.PareaTp2Hp50=Parea;
                     Cellrecord{fn}.chgTp2Hp50=chg;
                      Cellrecord{fn}.peakTp2Hp50=peak;
                          Cellrecord{fn}.MchgTp2Hp50=Mchg;
                      Cellrecord{fn}.MpeakTp2Hp50=Mpeak;
                      % Peak------------------------------------------
                    dx =Cellrecord{fn}.PatternSpacing(1); dy = dx;
                    xlim = [ceil(min(X))-dx-0.1 ceil(max(X))]; ylim = [ceil(min(Y))-dy-0.1 ceil(max(Y))];
                    
                    
                    xtick = xlim(1):dx:xlim(2); ytick = ylim(1):dy:ylim(2);
                    xi = xtick(1:end-1)+dx/2; yi = ytick(1:end-1)+dy/2;
                    avg = gridavg2(X,Y,Zpk,xtick,ytick,false,1,0);
                                       yy=nansum(avg,2);
                                       xx=nansum(avg);
                    [mxx ii]=max(xx);
                    jj=find(xx>mxx*Perc);
                    
                    if ~isempty(jj)
                        %                      ctext=abs(xi(jj(1))-xi(jj(end)));
                        ctextratio=abs(xi(jj(1))-xi(jj(end)))/abs(xlim(2)-xlim(1));
                        Dx=max(dx,abs(xi(jj(1))-xi(jj(end))));
                    else Dx=0;ctext=0;ctextratio=0;
                    end
                    [myy ii]=max(yy);
                    jj=find(yy>myy*Perc);
                    if ~isempty(jj)
                        Dy=max(dx,abs(yi(jj(1))-yi(jj(end))));
                    else Dy=0;
                    end
                  
                    D=[];Skpeak=[];Laydist=[];
                    for bi=1:length(Bnd)-1
                        
                            if bi==length(Bnd)-1
                             indBnd=find(xi>=Bnd(bi));
                        else
                            indBnd=find(xi>=Bnd(bi)&xi<Bnd(bi+1));
                        end
                        
                       
                        y1=nansum(avg(:,indBnd),2);
                        s1=skewness(y1);
                        Laydist=[Laydist y1];
                     
                        [mx ii]=max(y1);
                        jj=find(y1>mx*Perc);
                        
                        if ~isempty(jj)
                            D1=max(dx,abs(yi(jj(1))-yi(jj(end))));
                        else D1=0;
                        end
                        D=[D D1];
                        Skpeak=[Skpeak s1];
                    end
                   
                    DW=[Dx ctextratio  Dy Dy./abs(ylim(1)-ylim(2)) D];
                    DistWidthPeak=[DistWidthPeak;DW];
                    SkPeaklayer=[SkPeaklayer;Skpeak];
                    Cellrecord{fn}.LaydistPkTp2Hp50=Laydist;
                    Cellrecord{fn}.DistWidthPeakTp2Hp50=DW;
                    Cellrecord{fn}.SkPeaklayerTp2Hp50=Skpeak;
                   
                    
                     %charge----------------------------------------------------------------------------------
                    avg = gridavg2(X,Y,Zchg,xtick,ytick,false,1,0);
                    xx=nansum(avg);
                    yy=nansum(avg,2);
                    [mxx ii]=max(xx);
                    jj=find(xx>mxx*Perc);
                    
                    if ~isempty(jj)
                        %                      ctext=abs(xi(jj(1))-xi(jj(end)));
                        ctextratio=abs(xi(jj(1))-xi(jj(end)))/abs(xlim(2)-xlim(1));
                        Dx=max(dx,abs(xi(jj(1))-xi(jj(end))));
                    else Dx=0;ctext=0;ctextratio=0;
                    end
                    
                    [myy ii]=max(yy);
                    jj=find(yy>myy*Perc);
                    if ~isempty(jj)
                        Dy=max(dx,abs(yi(jj(1))-yi(jj(end))));
                    else Dy=0;
                    end
             
                    D=[];Skchg=[];Laydist=[];
                    for bi=1:length(Bnd)-1
                      
                            if bi==length(Bnd)-1
                             indBnd=find(xi>=Bnd(bi));
                        else
                            indBnd=find(xi>=Bnd(bi)&xi<Bnd(bi+1));
                        end
                        
                       
                        y1=nansum(avg(:,indBnd),2);
                        s1=skewness(y1);
                        Laydist=[Laydist y1];
                   
                        [mx ii]=max(y1);
                        jj=find(y1>mx*Perc);
                        
                        if ~isempty(jj)
                            D1=max(dx,abs(yi(jj(1))-yi(jj(end))));
                        else D1=0;
                        end
                        D=[D D1];
                        Skchg=[Skchg s1];
                    end
                    
                    DW=[Dx ctextratio  Dy Dy./abs(ylim(1)-ylim(2)) D];
                    DistWidthCharge=[DistWidthCharge;DW];
                  
                    Skchglayer=[Skchglayer;Skchg];
                    Cellrecord{fn}.LaydistChgTp2Hp50=Laydist;                    
                    Cellrecord{fn}.DistWidthChargeTp2Hp50=DW;
                    Cellrecord{fn}.SkchglayerTp2Hp50=Skchg;
                    
                    
                    %density---------------------------------------------------------------
                    avg = gridavg2(X,Y,Zflag,xtick,ytick,false,1,0);
                    xx=nansum(avg);
                    yy=nansum(avg,2);
                    [mxx ii]=max(xx);
                    jj=find(xx>mxx*Perc);
                    
                    if ~isempty(jj)
                        %                      ctext=abs(xi(jj(1))-xi(jj(end)));
                        ctextratio=abs(xi(jj(1))-xi(jj(end)))/abs(xlim(2)-xlim(1));
                        Dx=max(dx,abs(xi(jj(1))-xi(jj(end))));
                    else Dx=0;ctext=0;ctextratio=0;
                    end
                    [myy ii]=max(yy);
                    jj=find(yy>myy*Perc);
                    if ~isempty(jj)
                        Dy=max(dx,abs(yi(jj(1))-yi(jj(end))));
                    else Dy=0;
                    end
                    
                    D=[];Skdens=[];Laydist=[];
                    for bi=1:length(Bnd)-1
                        
                            if bi==length(Bnd)-1
                             indBnd=find(xi>=Bnd(bi));
                        else
                            indBnd=find(xi>=Bnd(bi)&xi<Bnd(bi+1));
                        end
                        
                      
                        
                        y1=nansum(avg(:,indBnd),2);
                        s1=skewness(y1);
                        Laydist=[Laydist y1];
                        
                        [mx ii]=max(y1);
                        jj=find(y1>mx*Perc);
                        
                        if ~isempty(jj)
                            D1=max(dx,abs(yi(jj(1))-yi(jj(end))));
                        else D1=0;
                        end
                        D=[D D1];
                        Skdens=[Skdens s1];
                    end
                    DW=[Dx ctextratio  Dy Dy./abs(ylim(1)-ylim(2)) D];
                    DistWidthDense=[DistWidthDense;DW];
                    Cellrecord{fn}.LaydistDensTp2Hp50=Laydist;                 
                    Skdenslayer=[Skdenslayer;Skdens];
                    Cellrecord{fn}.DistWidthdensTp2Hp50=DW;
                    Cellrecord{fn}.SkdenslayerTp2Hp50=Skdens;
                    
                end
            end
            
            if handles.avgTp==2&&abs(handles.avgHoldingPotential-10)<1&&isfield(Cellrecord{fn},'meanlatencyPtxInd10')&&~isempty(Cellrecord{fn}.meanlatencyPtxInd10)
                
                Soma= Cellrecord{fn}.SomaCoordinates;
                flipimg=Cellrecord{fn}.flipimg;
                pth=char(Cellrecord{fn}.Pth);
                stimcoordinates=Cellrecord{fn}.StimCoordinates;
                rotate_angle=Cellrecord{fn}.SpatialRotation;
                [stimCoordinates]=Matrixrotate_singlemap_GUI(Soma,stimcoordinates,flipimg,rotate_angle);
                X=stimCoordinates(1,:);
                Y=stimCoordinates(2,:);
                Zpk=zeros(size(X)).*NaN;
                Zchg=zeros(size(X)).*NaN;
                Zlat=zeros(size(X)).*NaN;
                indb=handles.indb;
                if abs(Bnd(indb)-Bnd(1))>=handles.mindist&&abs(Bnd(indb)-Bnd(1))<=handles.maxdist
                     N=N+1;
                xaxis=[xaxis,X];
                yaxis=[yaxis,Y];
                
                dX=[dX,Cellrecord{fn}.PatternSpacing(1)];
                Zflag=Cellrecord{fn}.meanflagPtxInd10;
                IndZ=find(Zflag>0);
                
                Zpk(IndZ)=Cellrecord{fn}.meanpeakPtxInd10(IndZ);
                
                Zchg(IndZ)=Cellrecord{fn}.meanareaPtxInd10(IndZ);
                
                Zlat(IndZ)=Cellrecord{fn}.meanlatencyPtxInd10(IndZ);
                   Pchg=[];
                    Ppeak=[];
                    chg=[];
                    peak=[];
                    A10=[];
                    Mchg=[];
                    Mpeak=[];
                    
                    for bi=1:length(Bnd)-1
                        
                            indBnd=find(X>=Bnd(bi)&X<Bnd(bi+1));
                        
                        
                        Pchg=[Pchg nansum(Zchg(indBnd))/nansum(Zchg)];
                        Ppeak=[Ppeak nansum(Zpk(indBnd))/nansum(Zpk)];
                        chg=[chg nansum(Zchg(indBnd))];
                        peak=[peak nansum(Zpk(indBnd))];
                         Mchg=[Mchg nanmean(Zchg(indBnd))];
                        Mpeak=[Mpeak nanmean(Zpk(indBnd))];
                        A10=[A10 nansum(Zpk(indBnd)>0).*Cellrecord{fn}.PatternSpacing(1).^2];
                    end
                    
                Pchglayer=[Pchglayer;Pchg];
                Ppeaklayer=[Ppeaklayer;Ppeak];
                Pk=[Pk,Zpk];
                
                
                Chg=[Chg, Zchg];
                
                Lat=[Lat,Zlat];
                Den=[Den Zflag'];
                Cellrecord{fn}.PchgTp2Hp10=Pchg;
                      Cellrecord{fn}.PpeakTp2Hp10=Ppeak;
                      Cellrecord{fn}.chgTp2Hp10=chg;
                      Cellrecord{fn}.peakTp2Hp10=peak;
                       Cellrecord{fn}.AreaTp2Hp10=A10;
                        Cellrecord{fn}.MchgTp2Hp10=Mchg;
                      Cellrecord{fn}.MpeakTp2Hp10=Mpeak;
                   % Peak------------------------------------------
                    dx =Cellrecord{fn}.PatternSpacing(1); dy = dx;
                    xlim = [ceil(min(X))-dx-0.1 ceil(max(X))]; ylim = [ceil(min(Y))-dy-0.1 ceil(max(Y))];
                    
                    
                    xtick = xlim(1):dx:xlim(2); ytick = ylim(1):dy:ylim(2);
                    xi = xtick(1:end-1)+dx/2; yi = ytick(1:end-1)+dy/2;
                    avg = gridavg2(X,Y,Zpk,xtick,ytick,false,1,0);
                                       yy=nansum(avg,2);
                                       xx=nansum(avg);
                    [mxx ii]=max(xx);
                    jj=find(xx>mxx*Perc);
                    
                    if ~isempty(jj)
                        %                      ctext=abs(xi(jj(1))-xi(jj(end)));
                        ctextratio=abs(xi(jj(1))-xi(jj(end)))/abs(xlim(2)-xlim(1));
                        Dx=max(dx,abs(xi(jj(1))-xi(jj(end))));
                    else Dx=0;ctext=0;ctextratio=0;
                    end
                    [myy ii]=max(yy);
                    jj=find(yy>myy*Perc);
                    if ~isempty(jj)
                        Dy=max(dx,abs(yi(jj(1))-yi(jj(end))));
                    else Dy=0;
                    end
                  
                    D=[];Skpeak=[];Laydist=[];
                    for bi=1:length(Bnd)-1
                     
                            if bi==length(Bnd)-1
                             indBnd=find(xi>=Bnd(bi));
                        else
                            indBnd=find(xi>=Bnd(bi)&xi<Bnd(bi+1));
                        end
                        
                       
                        
                        y1=nansum(avg(:,indBnd),2);
                        s1=skewness(y1);
                        Laydist=[Laydist y1];
                     
                        [mx ii]=max(y1);
                        jj=find(y1>mx*Perc);
                        
                        if ~isempty(jj)
                            D1=max(dx,abs(yi(jj(1))-yi(jj(end))));
                        else D1=0;
                        end
                        D=[D D1];
                        Skpeak=[Skpeak s1];
                    end
                   
                    DW=[Dx ctextratio  Dy Dy./abs(ylim(1)-ylim(2)) D];
                    DistWidthPeak=[DistWidthPeak;DW];
                    SkPeaklayer=[SkPeaklayer;Skpeak];
                    Cellrecord{fn}.LaydistPkTp2Hp10=Laydist;
                    Cellrecord{fn}.DistWidthPeakTp2Hp10=DW;
                    Cellrecord{fn}.SkPeaklayerTp2Hp10=Skpeak;
                    
                    
                    
                     %charge----------------------------------------------------------------------------------
                    avg = gridavg2(X,Y,Zchg,xtick,ytick,false,1,0);
                    xx=nansum(avg);
                    yy=nansum(avg,2);
                    [mxx ii]=max(xx);
                    jj=find(xx>mxx*Perc);
                    
                    if ~isempty(jj)
                        %                      ctext=abs(xi(jj(1))-xi(jj(end)));
                        ctextratio=abs(xi(jj(1))-xi(jj(end)))/abs(xlim(2)-xlim(1));
                        Dx=max(dx,abs(xi(jj(1))-xi(jj(end))));
                    else Dx=0;ctext=0;ctextratio=0;
                    end
                    
                    [myy ii]=max(yy);
                    jj=find(yy>myy*Perc);
                    if ~isempty(jj)
                        Dy=max(dx,abs(yi(jj(1))-yi(jj(end))));
                    else Dy=0;
                    end
             
                    D=[];Skchg=[];Laydist=[];
                    for bi=1:length(Bnd)-1
                        
                             if bi==length(Bnd)-1
                             indBnd=find(xi>=Bnd(bi));
                        else
                            indBnd=find(xi>=Bnd(bi)&xi<Bnd(bi+1));
                        end
                        
                      
                        
                        y1=nansum(avg(:,indBnd),2);
                        s1=skewness(y1);
                        Laydist=[Laydist y1];
                   
                        [mx ii]=max(y1);
                        jj=find(y1>mx*Perc);
                        
                        if ~isempty(jj)
                            D1=max(dx,abs(yi(jj(1))-yi(jj(end))));
                        else D1=0;
                        end
                        D=[D D1];
                        Skchg=[Skchg s1];
                    end
                    
                    DW=[Dx ctextratio  Dy Dy./abs(ylim(1)-ylim(2)) D];
                    DistWidthCharge=[DistWidthCharge;DW];
                  
                    Skchglayer=[Skchglayer;Skchg];
                    Cellrecord{fn}.LaydistChgTp2Hp10=Laydist;                    
                    Cellrecord{fn}.DistWidthChargeTp2Hp10=DW;
                    Cellrecord{fn}.SkchglayerTp2Hp10=Skchg;
                    
                    
                    %density---------------------------------------------------------------
                    avg = gridavg2(X,Y,Zflag,xtick,ytick,false,1,0);
                    xx=nansum(avg);
                    yy=nansum(avg,2);
                    [mxx ii]=max(xx);
                    jj=find(xx>mxx*Perc);
                    
                    if ~isempty(jj)
                        %                      ctext=abs(xi(jj(1))-xi(jj(end)));
                        ctextratio=abs(xi(jj(1))-xi(jj(end)))/abs(xlim(2)-xlim(1));
                        Dx=max(dx,abs(xi(jj(1))-xi(jj(end))));
                    else Dx=0;ctext=0;ctextratio=0;
                    end
                    [myy ii]=max(yy);
                    jj=find(yy>myy*Perc);
                    if ~isempty(jj)
                        Dy=max(dx,abs(yi(jj(1))-yi(jj(end))));
                    else Dy=0;
                    end
                    
                    D=[];Skdens=[];Laydist=[];
                    for bi=1:length(Bnd)-1
                      
                             if bi==length(Bnd)-1
                             indBnd=find(xi>=Bnd(bi));
                        else
                            indBnd=find(xi>=Bnd(bi)&xi<Bnd(bi+1));
                        end
                        
                        
                        y1=nansum(avg(:,indBnd),2);
                        s1=skewness(y1);
                        Laydist=[Laydist y1];
                        
                        [mx ii]=max(y1);
                        jj=find(y1>mx*Perc);
                        
                        if ~isempty(jj)
                            D1=max(dx,abs(yi(jj(1))-yi(jj(end))));
                        else D1=0;
                        end
                        D=[D D1];
                        Skdens=[Skdens s1];
                    end
                    DW=[Dx ctextratio  Dy Dy./abs(ylim(1)-ylim(2)) D];
                    DistWidthCharge=[DistWidthCharge;DW];
                    Cellrecord{fn}.LaydistChgTp2Hp10=Laydist;                 
                    Skdenslayer=[Skdenslayer;Skdens];
                    Cellrecord{fn}.DistWidthdensTp2Hp10=DW;
                    Cellrecord{fn}.SkdenslayerTp2Hp10=Skdens;
                    
                end
            end
            
            
            
            
            if handles.avgTp==3&&abs(handles.avgHoldingPotential+60)<1&&isfield(Cellrecord{fn},'meanlatencyapvInd60')&&~isempty(Cellrecord{fn}.meanlatencyapvInd60)
                
                Soma= Cellrecord{fn}.SomaCoordinates;
                flipimg=Cellrecord{fn}.flipimg;
                pth=char(Cellrecord{fn}.Pth);
                stimcoordinates=Cellrecord{fn}.StimCoordinates;
                rotate_angle=Cellrecord{fn}.SpatialRotation;
                indb=handles.indb;
                [stimCoordinates]=Matrixrotate_singlemap_GUI(Soma,stimcoordinates,flipimg,rotate_angle);
                X=stimCoordinates(1,:);
                Y=stimCoordinates(2,:);
                Zpk=zeros(size(X)).*NaN;
                Zchg=zeros(size(X)).*NaN;
                Zlat=zeros(size(X)).*NaN;
                if abs(Bnd(indb)-Bnd(1))>=handles.mindist&&abs(Bnd(indb)-Bnd(1))<=handles.maxdist
                     N=N+1;
                xaxis=[xaxis,X];
                yaxis=[yaxis,Y];
                
                dX=[dX,Cellrecord{fn}.PatternSpacing(1)];
                Zflag=Cellrecord{fn}.meanflagapvInd60;
                IndZ=find(Zflag>0);
                
                Zpk(IndZ)=Cellrecord{fn}.meanpeakapvInd60(IndZ);
                
                Zchg(IndZ)=Cellrecord{fn}.meanareaapvInd60(IndZ);
                
                Zlat(IndZ)=Cellrecord{fn}.meanlatencyapvInd60(IndZ);
                   Pchg=[];
                    Ppeak=[];
                    chg=[];
                    peak=[];
                    A60=[];
                     Mchg=[];
                    Mpeak=[];
                    for bi=1:length(Bnd)-1
                      
                            indBnd=find(X>=Bnd(bi)&X<Bnd(bi+1));
                     
                        
                        Pchg=[Pchg nansum(Zchg(indBnd))/nansum(Zchg)];
                        Ppeak=[Ppeak nansum(Zpk(indBnd))/nansum(Zpk)];
                        chg=[chg nansum(Zchg(indBnd))];
                        peak=[peak nansum(Zpk(indBnd))];
                         Mchg=[Mchg nanmean(Zchg(indBnd))];
                        Mpeak=[Mpeak nanmean(Zpk(indBnd))];
                        A60=[A60 nansum(Zpk(indBnd)>0).*Cellrecord{fn}.PatternSpacing(1).^2];
                    end
                    
                Pchglayer=[Pchglayer;Pchg];
                Ppeaklayer=[Ppeaklayer;Ppeak];
                Pk=[Pk,Zpk];
                
                
                Chg=[Chg, Zchg];
                
                Lat=[Lat,Zlat];
                Den=[Den Zflag'];
                Cellrecord{fn}.PchgTp3Hp60=Pchg;
                      Cellrecord{fn}.PpeakTp3Hp60=Ppeak;
                      Cellrecord{fn}.chgTp3Hp60=chg;
                      Cellrecord{fn}.peakTp3Hp60=peak;
                       Cellrecord{fn}.AreaTp3Hp60=A60;
                       Cellrecord{fn}.MchgTp3Hp60=Mchg;
                      Cellrecord{fn}.MpeakTp3Hp60=Mpeak;
                       
                     % Peak------------------------------------------
                    dx =Cellrecord{fn}.PatternSpacing(1); dy = dx;
                    xlim = [ceil(min(X))-dx-0.1 ceil(max(X))]; ylim = [ceil(min(Y))-dy-0.1 ceil(max(Y))];
                    
                    
                    xtick = xlim(1):dx:xlim(2); ytick = ylim(1):dy:ylim(2);
                    xi = xtick(1:end-1)+dx/2; yi = ytick(1:end-1)+dy/2;
                    avg = gridavg2(X,Y,Zpk,xtick,ytick,false,1,0);
                                       yy=nansum(avg,2);
                                       xx=nansum(avg);
                    [mxx ii]=max(xx);
                    jj=find(xx>mxx*Perc);
                    
                    if ~isempty(jj)
                        %                      ctext=abs(xi(jj(1))-xi(jj(end)));
                        ctextratio=abs(xi(jj(1))-xi(jj(end)))/abs(xlim(2)-xlim(1));
                        Dx=max(dx,abs(xi(jj(1))-xi(jj(end))));
                    else Dx=0;ctext=0;ctextratio=0;
                    end
                    [myy ii]=max(yy);
                    jj=find(yy>myy*Perc);
                    if ~isempty(jj)
                        Dy=max(dx,abs(yi(jj(1))-yi(jj(end))));
                    else Dy=0;
                    end
                  
                    D=[];Skpeak=[];Laydist=[];
                    for bi=1:length(Bnd)-1
                       
                             if bi==length(Bnd)-1
                             indBnd=find(xi>=Bnd(bi));
                        else
                            indBnd=find(xi>=Bnd(bi)&xi<Bnd(bi+1));
                        end
                        
                        
                        y1=nansum(avg(:,indBnd),2);
                        s1=skewness(y1);
                        Laydist=[Laydist y1];
                     
                        [mx ii]=max(y1);
                        jj=find(y1>mx*Perc);
                        
                        if ~isempty(jj)
                            D1=max(dx,abs(yi(jj(1))-yi(jj(end))));
                        else D1=0;
                        end
                        D=[D D1];
                        Skpeak=[Skpeak s1];
                    end
                   
                    DW=[Dx ctextratio  Dy Dy./abs(ylim(1)-ylim(2)) D];
                    DistWidthPeak=[DistWidthPeak;DW];
                    SkPeaklayer=[SkPeaklayer;Skpeak];
                    Cellrecord{fn}.LaydistPkTp3Hp60=Laydist;
                    Cellrecord{fn}.DistWidthPeakTp3Hp60=DW;
                    Cellrecord{fn}.SkPeaklayerTp3Hp60=Skpeak;
                    
                     ZflagD=Cellrecord{fn}.meanflagDapvInd60;
                    IndZD=find(ZflagD>0); 
                    Cellrecord{fn}.DirectAreaTp3Hp60=length(IndZD).*Cellrecord{fn}.PatternSpacing(1).^2;
                    DirArea=[DirArea;length(IndZD).*Cellrecord{fn}.PatternSpacing(1).^2];
                     %charge----------------------------------------------------------------------------------
                    avg = gridavg2(X,Y,Zchg,xtick,ytick,false,1,0);
                    xx=nansum(avg);
                    yy=nansum(avg,2);
                    [mxx ii]=max(xx);
                    jj=find(xx>mxx*Perc);
                    
                    if ~isempty(jj)
                        %                      ctext=abs(xi(jj(1))-xi(jj(end)));
                        ctextratio=abs(xi(jj(1))-xi(jj(end)))/abs(xlim(2)-xlim(1));
                        Dx=max(dx,abs(xi(jj(1))-xi(jj(end))));
                    else Dx=0;ctext=0;ctextratio=0;
                    end
                    
                    [myy ii]=max(yy);
                    jj=find(yy>myy*Perc);
                    if ~isempty(jj)
                        Dy=max(dx,abs(yi(jj(1))-yi(jj(end))));
                    else Dy=0;
                    end
             
                    D=[];Skchg=[];Laydist=[];
                    for bi=1:length(Bnd)-1
                       
                            if bi==length(Bnd)-1
                             indBnd=find(xi>=Bnd(bi));
                        else
                            indBnd=find(xi>=Bnd(bi)&xi<Bnd(bi+1));
                        end
                        
                        
                        
                        y1=nansum(avg(:,indBnd),2);
                        s1=skewness(y1);
                        Laydist=[Laydist y1];
                   
                        [mx ii]=max(y1);
                        jj=find(y1>mx*Perc);
                        
                        if ~isempty(jj)
                            D1=max(dx,abs(yi(jj(1))-yi(jj(end))));
                        else D1=0;
                        end
                        D=[D D1];
                        Skchg=[Skchg s1];
                    end
                    
                    DW=[Dx ctextratio  Dy Dy./abs(ylim(1)-ylim(2)) D];
                    DistWidthCharge=[DistWidthCharge;DW];
                  
                    Skchglayer=[Skchglayer;Skchg];
                    Cellrecord{fn}.LaydistChgTp3Hp60=Laydist;                    
                    Cellrecord{fn}.DistWidthChargeTp3Hp60=DW;
                    Cellrecord{fn}.SkchglayerTp3Hp60=Skchg;
                    
                    
                    %density---------------------------------------------------------------
                    avg = gridavg2(X,Y,Zflag,xtick,ytick,false,1,0);
                    xx=nansum(avg);
                    yy=nansum(avg,2);
                    [mxx ii]=max(xx);
                    jj=find(xx>mxx*Perc);
                    
                    if ~isempty(jj)
                        %                      ctext=abs(xi(jj(1))-xi(jj(end)));
                        ctextratio=abs(xi(jj(1))-xi(jj(end)))/abs(xlim(2)-xlim(1));
                        Dx=max(dx,abs(xi(jj(1))-xi(jj(end))));
                    else Dx=0;ctext=0;ctextratio=0;
                    end
                    [myy ii]=max(yy);
                    jj=find(yy>myy*Perc);
                    if ~isempty(jj)
                        Dy=max(dx,abs(yi(jj(1))-yi(jj(end))));
                    else Dy=0;
                    end
                    
                    D=[];Skdens=[];Laydist=[];
                    for bi=1:length(Bnd)-1
                        
                         if bi==length(Bnd)-1
                             indBnd=find(xi>=Bnd(bi));
                        else
                            indBnd=find(xi>=Bnd(bi)&xi<Bnd(bi+1));
                        end
                        
                        
                        
                        y1=nansum(avg(:,indBnd),2);
                        s1=skewness(y1);
                        Laydist=[Laydist y1];
                        
                        [mx ii]=max(y1);
                        jj=find(y1>mx*Perc);
                        
                        if ~isempty(jj)
                            D1=max(dx,abs(yi(jj(1))-yi(jj(end))));
                        else D1=0;
                        end
                        D=[D D1];
                        Skdens=[Skdens s1];
                    end
                    DW=[Dx ctextratio  Dy Dy./abs(ylim(1)-ylim(2)) D];
                    DistWidthCharge=[DistWidthCharge;DW];
                    Cellrecord{fn}.LaydistChgTp3Hp60=Laydist;                 
                    Skdenslayer=[Skdenslayer;Skdens];
                    Cellrecord{fn}.DistWidthdensTp3Hp60=DW;
                    Cellrecord{fn}.SkdenslayerTp3Hp60=Skdens;
                    
                end
            end
            
            
            
            
            if handles.avgTp==3&&abs(handles.avgHoldingPotential-50)<1&&isfield(Cellrecord{fn},'meanlatencyapvInd50')&&~isempty(Cellrecord{fn}.meanlatencyapvInd50)
              
                Soma= Cellrecord{fn}.SomaCoordinates;
                flipimg=Cellrecord{fn}.flipimg;
                pth=char(Cellrecord{fn}.Pth);
                stimcoordinates=Cellrecord{fn}.StimCoordinates;
                rotate_angle=Cellrecord{fn}.SpatialRotation;
                [stimCoordinates]=Matrixrotate_singlemap_GUI(Soma,stimcoordinates,flipimg,rotate_angle);
                X=stimCoordinates(1,:);
                Y=stimCoordinates(2,:);
                Zpk=zeros(size(X)).*NaN;
                Zchg=zeros(size(X)).*NaN;
                Zlat=zeros(size(X)).*NaN;
                indb=handles.indb;
                if abs(Bnd(indb)-Bnd(1))>=handles.mindist&&abs(Bnd(indb)-Bnd(1))<=handles.maxdist
                     N=N+1;
                xaxis=[xaxis,X];
                yaxis=[yaxis,Y];
             
                dX=[dX,Cellrecord{fn}.PatternSpacing(1)];
                Zflag=Cellrecord{fn}.meanflagapvInd50;
                IndZ=find(Zflag>0);
                
                Zpk(IndZ)=Cellrecord{fn}.meanpeakapvInd50(IndZ);
                
                Zchg(IndZ)=Cellrecord{fn}.meanareaapvInd50(IndZ);
                
                Zlat(IndZ)=Cellrecord{fn}.meanlatencyapvInd50(IndZ);
               Pchg=[];
                    Ppeak=[];
                    chg=[];
                    peak=[];
                    A50=[];
                    Mchg=[];
                    Mpeak=[];
                    for bi=1:length(Bnd)-1
                       
                            indBnd=find(X>=Bnd(bi)&X<Bnd(bi+1));
                     
                        
                        Pchg=[Pchg nansum(Zchg(indBnd))/nansum(Zchg)];
                        Ppeak=[Ppeak nansum(Zpk(indBnd))/nansum(Zpk)];
                        chg=[chg nansum(Zchg(indBnd))];
                        peak=[peak nansum(Zpk(indBnd))];
                          Mchg=[Mchg nanmean(Zchg(indBnd))];
                        Mpeak=[Mpeak nanmean(Zpk(indBnd))];
                        A50=[A50 nansum(Zpk(indBnd)>0).*Cellrecord{fn}.PatternSpacing(1).^2];
                    end
                Pchglayer=[Pchglayer;Pchg];
                Ppeaklayer=[Ppeaklayer;Ppeak];
                Pk=[Pk,Zpk];
                
                
                Chg=[Chg, Zchg];
                
                Lat=[Lat,Zlat];
                Den=[Den Zflag'];
                Cellrecord{fn}.PchgTp3Hp50=Pchg;
                      Cellrecord{fn}.PpeakTp3Hp50=Ppeak;
                      Cellrecord{fn}.chgTp3Hp50=chg;
                      Cellrecord{fn}.peakTp3Hp50=peak;
                       Cellrecord{fn}.AreaTp3Hp50=A50;
                       Cellrecord{fn}.MchgTp3Hp50=Mchg;
                      Cellrecord{fn}.MpeakTp3Hp50=Mpeak;
              % Peak------------------------------------------
                    dx =Cellrecord{fn}.PatternSpacing(1); dy = dx;
                    xlim = [ceil(min(X))-dx-0.1 ceil(max(X))]; ylim = [ceil(min(Y))-dy-0.1 ceil(max(Y))];
                    
                    
                    xtick = xlim(1):dx:xlim(2); ytick = ylim(1):dy:ylim(2);
                    xi = xtick(1:end-1)+dx/2; yi = ytick(1:end-1)+dy/2;
                    avg = gridavg2(X,Y,Zpk,xtick,ytick,false,1,0);
                                       yy=nansum(avg,2);
                                       xx=nansum(avg);
                    [mxx ii]=max(xx);
                    jj=find(xx>mxx*Perc);
                    
                    if ~isempty(jj)
                        %                      ctext=abs(xi(jj(1))-xi(jj(end)));
                        ctextratio=abs(xi(jj(1))-xi(jj(end)))/abs(xlim(2)-xlim(1));
                        Dx=max(dx,abs(xi(jj(1))-xi(jj(end))));
                    else Dx=0;ctext=0;ctextratio=0;
                    end
                    [myy ii]=max(yy);
                    jj=find(yy>myy*Perc);
                    if ~isempty(jj)
                        Dy=max(dx,abs(yi(jj(1))-yi(jj(end))));
                    else Dy=0;
                    end
                  
                    D=[];Skpeak=[];Laydist=[];
                    for bi=1:length(Bnd)-1
                        
                            if bi==length(Bnd)-1
                             indBnd=find(xi>=Bnd(bi));
                        else
                            indBnd=find(xi>=Bnd(bi)&xi<Bnd(bi+1));
                        end
                        
                     
                        
                        y1=nansum(avg(:,indBnd),2);
                        s1=skewness(y1);
                        Laydist=[Laydist y1];
                     
                        [mx ii]=max(y1);
                        jj=find(y1>mx*Perc);
                        
                        if ~isempty(jj)
                            D1=max(dx,abs(yi(jj(1))-yi(jj(end))));
                        else D1=0;
                        end
                        D=[D D1];
                        Skpeak=[Skpeak s1];
                    end
                   
                    DW=[Dx ctextratio  Dy Dy./abs(ylim(1)-ylim(2)) D];
                    DistWidthPeak=[DistWidthPeak;DW];
                    SkPeaklayer=[SkPeaklayer;Skpeak];
                    Cellrecord{fn}.LaydistPkTp3Hp50=Laydist;
                    Cellrecord{fn}.DistWidthPeakTp3Hp50=DW;
                    Cellrecord{fn}.SkPeaklayerTp3Hp50=Skpeak;
                    
                    
                    
                     %charge----------------------------------------------------------------------------------
                    avg = gridavg2(X,Y,Zchg,xtick,ytick,false,1,0);
                    xx=nansum(avg);
                    yy=nansum(avg,2);
                    [mxx ii]=max(xx);
                    jj=find(xx>mxx*Perc);
                    
                    if ~isempty(jj)
                        %                      ctext=abs(xi(jj(1))-xi(jj(end)));
                        ctextratio=abs(xi(jj(1))-xi(jj(end)))/abs(xlim(2)-xlim(1));
                        Dx=max(dx,abs(xi(jj(1))-xi(jj(end))));
                    else Dx=0;ctext=0;ctextratio=0;
                    end
                    
                    [myy ii]=max(yy);
                    jj=find(yy>myy*Perc);
                    if ~isempty(jj)
                        Dy=max(dx,abs(yi(jj(1))-yi(jj(end))));
                    else Dy=0;
                    end
             
                    D=[];Skchg=[];Laydist=[];
                    for bi=1:length(Bnd)-1
                       
                            if bi==length(Bnd)-1
                             indBnd=find(xi>=Bnd(bi));
                        else
                            indBnd=find(xi>=Bnd(bi)&xi<Bnd(bi+1));
                        end
                        
                       
                        
                        y1=nansum(avg(:,indBnd),2);
                        s1=skewness(y1);
                        Laydist=[Laydist y1];
                   
                        [mx ii]=max(y1);
                        jj=find(y1>mx*Perc);
                        
                        if ~isempty(jj)
                            D1=max(dx,abs(yi(jj(1))-yi(jj(end))));
                        else D1=0;
                        end
                        D=[D D1];
                        Skchg=[Skchg s1];
                    end
                    
                    DW=[Dx ctextratio  Dy Dy./abs(ylim(1)-ylim(2)) D];
                    DistWidthCharge=[DistWidthCharge;DW];
                  
                    Skchglayer=[Skchglayer;Skchg];
                    Cellrecord{fn}.LaydistChgTp3Hp50=Laydist;                    
                    Cellrecord{fn}.DistWidthChargeTp3Hp50=DW;
                    Cellrecord{fn}.SkchglayerTp3Hp50=Skchg;
                    
                    
                    %density---------------------------------------------------------------
                    avg = gridavg2(X,Y,Zflag,xtick,ytick,false,1,0);
                    xx=nansum(avg);
                    yy=nansum(avg,2);
                    [mxx ii]=max(xx);
                    jj=find(xx>mxx*Perc);
                    
                    if ~isempty(jj)
                        %                      ctext=abs(xi(jj(1))-xi(jj(end)));
                        ctextratio=abs(xi(jj(1))-xi(jj(end)))/abs(xlim(2)-xlim(1));
                        Dx=max(dx,abs(xi(jj(1))-xi(jj(end))));
                    else Dx=0;ctext=0;ctextratio=0;
                    end
                    [myy ii]=max(yy);
                    jj=find(yy>myy*Perc);
                    if ~isempty(jj)
                        Dy=max(dx,abs(yi(jj(1))-yi(jj(end))));
                    else Dy=0;
                    end
                    
                    D=[];Skdens=[];Laydist=[];
                    for bi=1:length(Bnd)-1
                       
                           if bi==length(Bnd)-1
                             indBnd=find(xi>=Bnd(bi));
                        else
                            indBnd=find(xi>=Bnd(bi)&xi<Bnd(bi+1));
                        end
                        
                       
                        
                        y1=nansum(avg(:,indBnd),2);
                        s1=skewness(y1);
                        Laydist=[Laydist y1];
                        
                        [mx ii]=max(y1);
                        jj=find(y1>mx*Perc);
                        
                        if ~isempty(jj)
                            D1=max(dx,abs(yi(jj(1))-yi(jj(end))));
                        else D1=0;
                        end
                        D=[D D1];
                        Skdens=[Skdens s1];
                    end
                    DW=[Dx ctextratio  Dy Dy./abs(ylim(1)-ylim(2)) D];
                    DistWidthCharge=[DistWidthCharge;DW];
                    Cellrecord{fn}.LaydistChgTp3Hp50=Laydist;                 
                    Skdenslayer=[Skdenslayer;Skdens];
                    Cellrecord{fn}.DistWidthdensTp3Hp50=DW;
                    Cellrecord{fn}.SkdenslayerTp3Hp50=Skdens;
                    
                end
            end
            
            
            
            if handles.avgTp==3&&abs(handles.avgHoldingPotential-10)<1&&~isempty(Cellrecord{fn}.meanlatencyapvInd10)
               
                Soma= Cellrecord{fn}.SomaCoordinates;
                flipimg=Cellrecord{fn}.flipimg;
                pth=char(Cellrecord{fn}.Pth);
                stimcoordinates=Cellrecord{fn}.StimCoordinates;
                rotate_angle=Cellrecord{fn}.SpatialRotation;
                [stimCoordinates]=Matrixrotate_singlemap_GUI(Soma,stimcoordinates,flipimg,rotate_angle);
                X=stimCoordinates(1,:);
                Y=stimCoordinates(2,:);
                Zpk=zeros(size(X)).*NaN;
                Zchg=zeros(size(X)).*NaN;
                Zlat=zeros(size(X)).*NaN;
                indb=handles.indb;
                if abs(Bnd(indb)-Bnd(1))>=handles.mindist&&abs(Bnd(indb)-Bnd(1))<=handles.maxdist
                     N=N+1;
                xaxis=[xaxis,X];
                yaxis=[yaxis,Y];
               
                dX=[dX,Cellrecord{fn}.PatternSpacing(1)];
                Zflag=Cellrecord{fn}.meanflagapvInd10;
                IndZ=find(Zflag>0);
                
                Zpk(IndZ)=Cellrecord{fn}.meanpeakapvInd10(IndZ);
                
                Zchg(IndZ)=Cellrecord{fn}.meanareaapvInd10(IndZ);
                
                Zlat(IndZ)=Cellrecord{fn}.meanlatencyapvInd10(IndZ);
                Pchg=[];
                    Ppeak=[];
                    chg=[];
                    peak=[];
                    A10=[];
                     Mchg=[];
                    Mpeak=[];
                    
                    for bi=1:length(Bnd)-1
                        
                            indBnd=find(X>=Bnd(bi)&X<Bnd(bi+1));
                       
                        
                        Pchg=[Pchg nansum(Zchg(indBnd))/nansum(Zchg)];
                        Ppeak=[Ppeak nansum(Zpk(indBnd))/nansum(Zpk)];
                        chg=[chg nansum(Zchg(indBnd))];
                        peak=[peak nansum(Zpk(indBnd))];
                        Mchg=[Mchg nanmean(Zchg(indBnd))];
                        Mpeak=[Mpeak nanmean(Zpk(indBnd))];
                        A10=[A10 nansum(Zpk(indBnd)>0).*Cellrecord{fn}.PatternSpacing(1).^2];
                    end
                Pchglayer=[Pchglayer;Pchg];
                Ppeaklayer=[Ppeaklayer;Ppeak];
                Pk=[Pk,Zpk];
                
                
                Chg=[Chg, Zchg];
                
                Lat=[Lat,Zlat];
                Den=[Den Zflag'];
                Cellrecord{fn}.PchgTp3Hp10=Pchg;
                      Cellrecord{fn}.PpeakTp3Hp10=Ppeak;
                      Cellrecord{fn}.chgTp3Hp10=chg;
                      Cellrecord{fn}.peakTp3Hp10=peak;
                       Cellrecord{fn}.AreaTp3Hp10=A10;
                       Cellrecord{fn}.MchgTp3Hp10=Mchg;
                      Cellrecord{fn}.MpeakTp3Hp10=Mpeak;
                   % Peak------------------------------------------
                    dx =Cellrecord{fn}.PatternSpacing(1); dy = dx;
                    xlim = [ceil(min(X))-dx-0.1 ceil(max(X))]; ylim = [ceil(min(Y))-dy-0.1 ceil(max(Y))];
                    
                    
                    xtick = xlim(1):dx:xlim(2); ytick = ylim(1):dy:ylim(2);
                    xi = xtick(1:end-1)+dx/2; yi = ytick(1:end-1)+dy/2;
                    avg = gridavg2(X,Y,Zpk,xtick,ytick,false,1,0);
                                       yy=nansum(avg,2);
                                       xx=nansum(avg);
                    [mxx ii]=max(xx);
                    jj=find(xx>mxx*Perc);
                    
                    if ~isempty(jj)
                        %                      ctext=abs(xi(jj(1))-xi(jj(end)));
                        ctextratio=abs(xi(jj(1))-xi(jj(end)))/abs(xlim(2)-xlim(1));
                        Dx=max(dx,abs(xi(jj(1))-xi(jj(end))));
                    else Dx=0;ctext=0;ctextratio=0;
                    end
                    [myy ii]=max(yy);
                    jj=find(yy>myy*Perc);
                    if ~isempty(jj)
                        Dy=max(dx,abs(yi(jj(1))-yi(jj(end))));
                    else Dy=0;
                    end
                  
                    D=[];Skpeak=[];Laydist=[];
                    for bi=1:length(Bnd)-1
                       
                        if bi==length(Bnd)-1
                             indBnd=find(xi>=Bnd(bi));
                        else
                            indBnd=find(xi>=Bnd(bi)&xi<Bnd(bi+1));
                        end
                        
                      
                        y1=nansum(avg(:,indBnd),2);
                        s1=skewness(y1);
                        Laydist=[Laydist y1];
                     
                        [mx ii]=max(y1);
                        jj=find(y1>mx*Perc);
                        
                        if ~isempty(jj)
                            D1=max(dx,abs(yi(jj(1))-yi(jj(end))));
                        else D1=0;
                        end
                        D=[D D1];
                        Skpeak=[Skpeak s1];
                    end
                   
                    DW=[Dx ctextratio  Dy Dy./abs(ylim(1)-ylim(2)) D];
                    DistWidthPeak=[DistWidthPeak;DW];
                    SkPeaklayer=[SkPeaklayer;Skpeak];
                    Cellrecord{fn}.LaydistPkTp3Hp10=Laydist;
                    Cellrecord{fn}.DistWidthPeakTp3Hp10=DW;
                    Cellrecord{fn}.SkPeaklayerTp3Hp10=Skpeak;
                    
                    
                    
                     %charge----------------------------------------------------------------------------------
                    avg = gridavg2(X,Y,Zchg,xtick,ytick,false,1,0);
                    xx=nansum(avg);
                    yy=nansum(avg,2);
                    [mxx ii]=max(xx);
                    jj=find(xx>mxx*Perc);
                    
                    if ~isempty(jj)
                        %                      ctext=abs(xi(jj(1))-xi(jj(end)));
                        ctextratio=abs(xi(jj(1))-xi(jj(end)))/abs(xlim(2)-xlim(1));
                        Dx=max(dx,abs(xi(jj(1))-xi(jj(end))));
                    else Dx=0;ctext=0;ctextratio=0;
                    end
                    
                    [myy ii]=max(yy);
                    jj=find(yy>myy*Perc);
                    if ~isempty(jj)
                        Dy=max(dx,abs(yi(jj(1))-yi(jj(end))));
                    else Dy=0;
                    end
             
                    D=[];Skchg=[];Laydist=[];
                    for bi=1:length(Bnd)-1
                        
                           if bi==length(Bnd)-1
                             indBnd=find(xi>=Bnd(bi));
                        else
                            indBnd=find(xi>=Bnd(bi)&xi<Bnd(bi+1));
                        end
                        
                      
                        
                        y1=nansum(avg(:,indBnd),2);
                        s1=skewness(y1);
                        Laydist=[Laydist y1];
                   
                        [mx ii]=max(y1);
                        jj=find(y1>mx*Perc);
                        
                        if ~isempty(jj)
                            D1=max(dx,abs(yi(jj(1))-yi(jj(end))));
                        else D1=0;
                        end
                        D=[D D1];
                        Skchg=[Skchg s1];
                    end
                    
                    DW=[Dx ctextratio  Dy Dy./abs(ylim(1)-ylim(2)) D];
                    DistWidthCharge=[DistWidthCharge;DW];
                  
                    Skchglayer=[Skchglayer;Skchg];
                    Cellrecord{fn}.LaydistChgTp3Hp10=Laydist;                    
                    Cellrecord{fn}.DistWidthChargeTp3Hp10=DW;
                    Cellrecord{fn}.SkchglayerTp3Hp10=Skchg;
                    
                    
                    %density---------------------------------------------------------------
                    avg = gridavg2(X,Y,Zflag,xtick,ytick,false,1,0);
                    xx=nansum(avg);
                    yy=nansum(avg,2);
                    [mxx ii]=max(xx);
                    jj=find(xx>mxx*Perc);
                    
                    if ~isempty(jj)
                        %                      ctext=abs(xi(jj(1))-xi(jj(end)));
                        ctextratio=abs(xi(jj(1))-xi(jj(end)))/abs(xlim(2)-xlim(1));
                        Dx=max(dx,abs(xi(jj(1))-xi(jj(end))));
                    else Dx=0;ctext=0;ctextratio=0;
                    end
                    [myy ii]=max(yy);
                    jj=find(yy>myy*Perc);
                    if ~isempty(jj)
                        Dy=max(dx,abs(yi(jj(1))-yi(jj(end))));
                    else Dy=0;
                    end
                    
                    D=[];Skdens=[];Laydist=[];
                    for bi=1:length(Bnd)-1
                 
                           if bi==length(Bnd)-1
                             indBnd=find(xi>=Bnd(bi));
                        else
                            indBnd=find(xi>=Bnd(bi)&xi<Bnd(bi+1));
                        end
                        
                      
                        y1=nansum(avg(:,indBnd),2);
                        s1=skewness(y1);
                        Laydist=[Laydist y1];
                        
                        [mx ii]=max(y1);
                        jj=find(y1>mx*Perc);
                        
                        if ~isempty(jj)
                            D1=max(dx,abs(yi(jj(1))-yi(jj(end))));
                        else D1=0;
                        end
                        D=[D D1];
                        Skdens=[Skdens s1];
                    end
                    DW=[Dx ctextratio  Dy Dy./abs(ylim(1)-ylim(2)) D];
                    DistWidthCharge=[DistWidthCharge;DW];
                    Cellrecord{fn}.LaydistChgTp3Hp10=Laydist;                 
                    Skdenslayer=[Skdenslayer;Skdens];
                    Cellrecord{fn}.DistWidthdensTp3Hp10=DW;
                    Cellrecord{fn}.SkdenslayerTp3Hp10=Skdens;
                    
                end
            end
            
            
            if handles.avgTp==4&&abs(handles.avgHoldingPotential+60)<1&&isfield(Cellrecord{fn},'meanlatencyttxInd60')&&~isempty(Cellrecord{fn}.meanlatencyttxInd60)
               
                Soma= Cellrecord{fn}.SomaCoordinates;
                flipimg=Cellrecord{fn}.flipimg;
                pth=char(Cellrecord{fn}.Pth);
                stimcoordinates=Cellrecord{fn}.StimCoordinates;
                rotate_angle=Cellrecord{fn}.SpatialRotation;
                [stimCoordinates]=Matrixrotate_singlemap_GUI(Soma,stimcoordinates,flipimg,rotate_angle);
                X=stimCoordinates(1,:);
                Y=stimCoordinates(2,:);
                Zpk=zeros(size(X)).*NaN;
                Zchg=zeros(size(X)).*NaN;
                Zlat=zeros(size(X)).*NaN;
                if abs(Bnd(indb)-Bnd(1))>=handles.mindist&&abs(Bnd(indb)-Bnd(1))<=handles.maxdist
                     N=N+1;
                xaxis=[xaxis,X];
                yaxis=[yaxis,Y];
                
                dX=[dX,Cellrecord{fn}.PatternSpacing(1)];
                Zflag=Cellrecord{fn}.meanflagttxInd60;
                IndZ=find(Zflag);
                
                Zpk(IndZ)=Cellrecord{fn}.meanpeakttxInd60(IndZ);
                
                Zchg(IndZ)=Cellrecord{fn}.meanareattxInd60(IndZ);
                
                Zlat(IndZ)=Cellrecord{fn}.meanlatencyttxInd60(IndZ);
                
                Pchg=[];
                    Ppeak=[];
                    chg=[];
                    peak=[];
                    A60=[];
                     Mchg=[];
                    Mpeak=[];
                    for bi=1:length(Bnd)-1
                     
                            indBnd=find(X>=Bnd(bi)&X<Bnd(bi+1));
                        
                        
                        Pchg=[Pchg nansum(Zchg(indBnd))/nansum(Zchg)];
                        Ppeak=[Ppeak nansum(Zpk(indBnd))/nansum(Zpk)];
                        chg=[chg nansum(Zchg(indBnd))];
                        peak=[peak nansum(Zpk(indBnd))];
                         Mchg=[Mchg nanmean(Zchg(indBnd))];
                        Mpeak=[Mpeak nanmean(Zpk(indBnd))];
                        A60=[A60 nansum(Zpk(indBnd)>0).*Cellrecord{fn}.PatternSpacing(1).^2];
                    end
                Pk=[Pk,Zpk];
                  Pchglayer=[Pchglayer;Pchg];
                Ppeaklayer=[Ppeaklayer;Ppeak];
                
                Chg=[Chg, Zchg];
                
                Lat=[Lat,Zlat];
                Den=[Den Zflag'];
                Cellrecord{fn}.PchgTp4Hp60=Pchg;
                Cellrecord{fn}.PpeakTp4Hp60=Ppeak;
                Cellrecord{fn}.chgTp4Hp60=chg;
                Cellrecord{fn}.peakTp4Hp60=peak;
                Cellrecord{fn}.AreaTp4Hp60=A60;
                
                Cellrecord{fn}.MchgTp4Hp60=Mchg;
                Cellrecord{fn}.MpeakTp4Hp60=Mpeak;
                ZflagD=Cellrecord{fn}.meanflagDttxInd60;
                IndZD=find(ZflagD>0);
                Cellrecord{fn}.DirectAreaTp4Hp60=length(IndZD).*Cellrecord{fn}.PatternSpacing(1).^2;
                DirArea=[DirArea;length(IndZD).*Cellrecord{fn}.PatternSpacing(1).^2];
                
                % Peak------------------------------------------
                dx =Cellrecord{fn}.PatternSpacing(1); dy = dx;
                    xlim = [ceil(min(X))-dx-0.1 ceil(max(X))]; ylim = [ceil(min(Y))-dy-0.1 ceil(max(Y))];
                    
                    
                    xtick = xlim(1):dx:xlim(2); ytick = ylim(1):dy:ylim(2);
                    xi = xtick(1:end-1)+dx/2; yi = ytick(1:end-1)+dy/2;
                    avg = gridavg2(X,Y,Zpk,xtick,ytick,false,1,0);
                                       yy=nansum(avg,2);
                                       xx=nansum(avg);
                    [mxx ii]=max(xx);
                    jj=find(xx>mxx*Perc);
                    
                    if ~isempty(jj)
                        %                      ctext=abs(xi(jj(1))-xi(jj(end)));
                        ctextratio=abs(xi(jj(1))-xi(jj(end)))/abs(xlim(2)-xlim(1));
                        Dx=max(dx,abs(xi(jj(1))-xi(jj(end))));
                    else Dx=0;ctext=0;ctextratio=0;
                    end
                    [myy ii]=max(yy);
                    jj=find(yy>myy*Perc);
                    if ~isempty(jj)
                        Dy=max(dx,abs(yi(jj(1))-yi(jj(end))));
                    else Dy=0;
                    end
                  
                    D=[];Skpeak=[];Laydist=[];
                    for bi=1:length(Bnd)-1
                      
                        if bi==length(Bnd)-1
                             indBnd=find(xi>=Bnd(bi));
                        else
                            indBnd=find(xi>=Bnd(bi)&xi<Bnd(bi+1));
                        end
                        
                       
                        
                        y1=nansum(avg(:,indBnd),2);
                        s1=skewness(y1);
                        Laydist=[Laydist y1];
                     
                        [mx ii]=max(y1);
                        jj=find(y1>mx*Perc);
                        
                        if ~isempty(jj)
                            D1=max(dx,abs(yi(jj(1))-yi(jj(end))));
                        else D1=0;
                        end
                        D=[D D1];
                        Skpeak=[Skpeak s1];
                    end
                   
                    DW=[Dx ctextratio  Dy Dy./abs(ylim(1)-ylim(2)) D];
                    DistWidthPeak=[DistWidthPeak;DW];
                    SkPeaklayer=[SkPeaklayer;Skpeak];
                    Cellrecord{fn}.LaydistPkTp4Hp60=Laydist;
                    Cellrecord{fn}.DistWidthPeakTp4Hp60=DW;
                    Cellrecord{fn}.SkPeaklayerTp4Hp60=Skpeak;
                    
                    
                    
                     %charge----------------------------------------------------------------------------------
                    avg = gridavg2(X,Y,Zchg,xtick,ytick,false,1,0);
                    xx=nansum(avg);
                    yy=nansum(avg,2);
                    [mxx ii]=max(xx);
                    jj=find(xx>mxx*Perc);
                    
                    if ~isempty(jj)
                        %                      ctext=abs(xi(jj(1))-xi(jj(end)));
                        ctextratio=abs(xi(jj(1))-xi(jj(end)))/abs(xlim(2)-xlim(1));
                        Dx=max(dx,abs(xi(jj(1))-xi(jj(end))));
                    else Dx=0;ctext=0;ctextratio=0;
                    end
                    
                    [myy ii]=max(yy);
                    jj=find(yy>myy*Perc);
                    if ~isempty(jj)
                        Dy=max(dx,abs(yi(jj(1))-yi(jj(end))));
                    else Dy=0;
                    end
             
                    D=[];Skchg=[];Laydist=[];
                    for bi=1:length(Bnd)-1
                    
                           if bi==length(Bnd)-1
                             indBnd=find(xi>=Bnd(bi));
                        else
                            indBnd=find(xi>=Bnd(bi)&xi<Bnd(bi+1));
                        end
                        
                        
                        
                        y1=nansum(avg(:,indBnd),2);
                        s1=skewness(y1);
                        Laydist=[Laydist y1];
                   
                        [mx ii]=max(y1);
                        jj=find(y1>mx*Perc);
                        
                        if ~isempty(jj)
                            D1=max(dx,abs(yi(jj(1))-yi(jj(end))));
                        else D1=0;
                        end
                        D=[D D1];
                        Skchg=[Skchg s1];
                    end
                    
                    DW=[Dx ctextratio  Dy Dy./abs(ylim(1)-ylim(2)) D];
                    DistWidthCharge=[DistWidthCharge;DW];
                  
                    Skchglayer=[Skchglayer;Skchg];
                    Cellrecord{fn}.LaydistChgTp4Hp60=Laydist;                    
                    Cellrecord{fn}.DistWidthChargeTp4Hp60=DW;
                    Cellrecord{fn}.SkchglayerTp4Hp60=Skchg;
                    
                    
                    %density---------------------------------------------------------------
                    avg = gridavg2(X,Y,Zflag,xtick,ytick,false,1,0);
                    xx=nansum(avg);
                    yy=nansum(avg,2);
                    [mxx ii]=max(xx);
                    jj=find(xx>mxx*Perc);
                    
                    if ~isempty(jj)
                        %                      ctext=abs(xi(jj(1))-xi(jj(end)));
                        ctextratio=abs(xi(jj(1))-xi(jj(end)))/abs(xlim(2)-xlim(1));
                        Dx=max(dx,abs(xi(jj(1))-xi(jj(end))));
                    else Dx=0;ctext=0;ctextratio=0;
                    end
                    [myy ii]=max(yy);
                    jj=find(yy>myy*Perc);
                    if ~isempty(jj)
                        Dy=max(dx,abs(yi(jj(1))-yi(jj(end))));
                    else Dy=0;
                    end
                    
                    D=[];Skdens=[];Laydist=[];
                    for bi=1:length(Bnd)-1
                        
                        if bi==length(Bnd)-1
                             indBnd=find(xi>=Bnd(bi));
                        else
                            indBnd=find(xi>=Bnd(bi)&xi<Bnd(bi+1));
                        end
                        
                        
                        y1=nansum(avg(:,indBnd),2);
                        s1=skewness(y1);
                        Laydist=[Laydist y1];
                        
                        [mx ii]=max(y1);
                        jj=find(y1>mx*Perc);
                        
                        if ~isempty(jj)
                            D1=max(dx,abs(yi(jj(1))-yi(jj(end))));
                        else D1=0;
                        end
                        D=[D D1];
                        Skdens=[Skdens s1];
                    end
                    DW=[Dx ctextratio  Dy Dy./abs(ylim(1)-ylim(2)) D];
                    DistWidthCharge=[DistWidthCharge;DW];
                    Cellrecord{fn}.LaydistChgTp4Hp60=Laydist;                 
                    Skdenslayer=[Skdenslayer;Skdens];
                    Cellrecord{fn}.DistWidthdensTp4Hp60=DW;
                    Cellrecord{fn}.SkdenslayerTp4Hp60=Skdens;
                    
                end
            end
            
            
            if handles.avgTp==4&&abs(handles.avgHoldingPotential-50)<1&&isfield(Cellrecord{fn},'meanlatencyttxInd50')&&~isempty(Cellrecord{fn}.meanlatencyttxInd50)
               
                Soma= Cellrecord{fn}.SomaCoordinates;
                flipimg=Cellrecord{fn}.flipimg;
                pth=char(Cellrecord{fn}.Pth);
                stimcoordinates=Cellrecord{fn}.StimCoordinates;
                rotate_angle=Cellrecord{fn}.SpatialRotation;
                [stimCoordinates]=Matrixrotate_singlemap_GUI(Soma,stimcoordinates,flipimg,rotate_angle);
                X=stimCoordinates(1,:);
                Y=stimCoordinates(2,:);
                Zpk=zeros(size(X)).*NaN;
                Zchg=zeros(size(X)).*NaN;
                Zlat=zeros(size(X)).*NaN;
                if abs(Bnd(indb)-Bnd(1))>=handles.mindist&&abs(Bnd(indb)-Bnd(1))<=handles.maxdist
                     N=N+1;
                xaxis=[xaxis,X];
                yaxis=[yaxis,Y];
                
                dX=[dX,Cellrecord{fn}.PatternSpacing(1)];
                Zflag=Cellrecord{fn}.meanflagttxInd50;
                IndZ=find(Zflag);
                
                Zpk(IndZ)=Cellrecord{fn}.meanpeakttxInd50(IndZ);
                
                Zchg(IndZ)=Cellrecord{fn}.meanareattxInd50(IndZ);
                
                Zlat(IndZ)=Cellrecord{fn}.meanlatencyttxInd50(IndZ);
               Pchg=[];
                    Ppeak=[];
                    chg=[];
                    peak=[];
                    A50=[];
                    Mchg=[];
                    Mpeak=[];
                    for bi=1:length(Bnd)-1
                       
                            indBnd=find(X>=Bnd(bi)&X<Bnd(bi+1));
                       
                        
                        Pchg=[Pchg nansum(Zchg(indBnd))/nansum(Zchg)];
                        Ppeak=[Ppeak nansum(Zpk(indBnd))/nansum(Zpk)];
                        chg=[chg nansum(Zchg(indBnd))];
                        peak=[peak nansum(Zpk(indBnd))];
                        Mchg=[Mchg nanmean(Zchg(indBnd))];
                        Mpeak=[Mpeak nanmean(Zpk(indBnd))];
                        A50=[A50 nansum(Zpk(indBnd)>0).*Cellrecord{fn}.PatternSpacing(1).^2];
                    end
                Pchglayer=[Pchglayer;Pchg];
                Ppeaklayer=[Ppeaklayer;Ppeak];
                Pk=[Pk,Zpk];
                
                
                Chg=[Chg, Zchg];
                
                Lat=[Lat,Zlat];
                Den=[Den Zflag'];
                Cellrecord{fn}.PchgTp4Hp50=Pchg;
                      Cellrecord{fn}.PpeakTp4Hp50=Ppeak;
                       Cellrecord{fn}.chgTp4Hp50=chg;
                      Cellrecord{fn}.peakTp4Hp50=peak;
                       Cellrecord{fn}.AreaTp4Hp50=A50;
                       Cellrecord{fn}.MchgTp4Hp50=Mchg;
                      Cellrecord{fn}.MpeakTp4Hp50=Mpeak;
              % Peak------------------------------------------
                    dx =Cellrecord{fn}.PatternSpacing(1); dy = dx;
                    xlim = [ceil(min(X))-dx-0.1 ceil(max(X))]; ylim = [ceil(min(Y))-dy-0.1 ceil(max(Y))];
                    
                    
                    xtick = xlim(1):dx:xlim(2); ytick = ylim(1):dy:ylim(2);
                    xi = xtick(1:end-1)+dx/2; yi = ytick(1:end-1)+dy/2;
                    avg = gridavg2(X,Y,Zpk,xtick,ytick,false,1,0);
                                       yy=nansum(avg,2);
                                       xx=nansum(avg);
                    [mxx ii]=max(xx);
                    jj=find(xx>mxx*Perc);
                    
                    if ~isempty(jj)
                        %                      ctext=abs(xi(jj(1))-xi(jj(end)));
                        ctextratio=abs(xi(jj(1))-xi(jj(end)))/abs(xlim(2)-xlim(1));
                        Dx=max(dx,abs(xi(jj(1))-xi(jj(end))));
                    else Dx=0;ctext=0;ctextratio=0;
                    end
                    [myy ii]=max(yy);
                    jj=find(yy>myy*Perc);
                    if ~isempty(jj)
                        Dy=max(dx,abs(yi(jj(1))-yi(jj(end))));
                    else Dy=0;
                    end
                  
                    D=[];Skpeak=[];Laydist=[];
                    for bi=1:length(Bnd)-1
                       
                            if bi==length(Bnd)-1
                             indBnd=find(xi>=Bnd(bi));
                        else
                            indBnd=find(xi>=Bnd(bi)&xi<Bnd(bi+1));
                        end
                        
                        
                        y1=nansum(avg(:,indBnd),2);
                        s1=skewness(y1);
                        Laydist=[Laydist y1];
                     
                        [mx ii]=max(y1);
                        jj=find(y1>mx*Perc);
                        
                        if ~isempty(jj)
                            D1=max(dx,abs(yi(jj(1))-yi(jj(end))));
                        else D1=0;
                        end
                        D=[D D1];
                        Skpeak=[Skpeak s1];
                    end
                   
                    DW=[Dx ctextratio  Dy Dy./abs(ylim(1)-ylim(2)) D];
                    DistWidthPeak=[DistWidthPeak;DW];
                    SkPeaklayer=[SkPeaklayer;Skpeak];
                    Cellrecord{fn}.LaydistPkTp4Hp50=Laydist;
                    Cellrecord{fn}.DistWidthPeakTp4Hp50=DW;
                    Cellrecord{fn}.SkPeaklayerTp4Hp50=Skpeak;
                    
                    
                    
                     %charge----------------------------------------------------------------------------------
                    avg = gridavg2(X,Y,Zchg,xtick,ytick,false,1,0);
                    xx=nansum(avg);
                    yy=nansum(avg,2);
                    [mxx ii]=max(xx);
                    jj=find(xx>mxx*Perc);
                    
                    if ~isempty(jj)
                        %                      ctext=abs(xi(jj(1))-xi(jj(end)));
                        ctextratio=abs(xi(jj(1))-xi(jj(end)))/abs(xlim(2)-xlim(1));
                        Dx=max(dx,abs(xi(jj(1))-xi(jj(end))));
                    else Dx=0;ctext=0;ctextratio=0;
                    end
                    
                    [myy ii]=max(yy);
                    jj=find(yy>myy*Perc);
                    if ~isempty(jj)
                        Dy=max(dx,abs(yi(jj(1))-yi(jj(end))));
                    else Dy=0;
                    end
             
                    D=[];Skchg=[];Laydist=[];
                    for bi=1:length(Bnd)-1
                      
                             if bi==length(Bnd)-1
                             indBnd=find(xi>=Bnd(bi));
                        else
                            indBnd=find(xi>=Bnd(bi)&xi<Bnd(bi+1));
                        end
                        
                   
                        
                        y1=nansum(avg(:,indBnd),2);
                        s1=skewness(y1);
                        Laydist=[Laydist y1];
                   
                        [mx ii]=max(y1);
                        jj=find(y1>mx*Perc);
                        
                        if ~isempty(jj)
                            D1=max(dx,abs(yi(jj(1))-yi(jj(end))));
                        else D1=0;
                        end
                        D=[D D1];
                        Skchg=[Skchg s1];
                    end
                    
                    DW=[Dx ctextratio  Dy Dy./abs(ylim(1)-ylim(2)) D];
                    DistWidthCharge=[DistWidthCharge;DW];
                  
                    Skchglayer=[Skchglayer;Skchg];
                    Cellrecord{fn}.LaydistChgTp4Hp50=Laydist;                    
                    Cellrecord{fn}.DistWidthChargeTp4Hp50=DW;
                    Cellrecord{fn}.SkchglayerTp4Hp50=Skchg;
                    
                    
                    %density---------------------------------------------------------------
                    avg = gridavg2(X,Y,Zflag,xtick,ytick,false,1,0);
                    xx=nansum(avg);
                    yy=nansum(avg,2);
                    [mxx ii]=max(xx);
                    jj=find(xx>mxx*Perc);
                    
                    if ~isempty(jj)
                        %                      ctext=abs(xi(jj(1))-xi(jj(end)));
                        ctextratio=abs(xi(jj(1))-xi(jj(end)))/abs(xlim(2)-xlim(1));
                        Dx=max(dx,abs(xi(jj(1))-xi(jj(end))));
                    else Dx=0;ctext=0;ctextratio=0;
                    end
                    [myy ii]=max(yy);
                    jj=find(yy>myy*Perc);
                    if ~isempty(jj)
                        Dy=max(dx,abs(yi(jj(1))-yi(jj(end))));
                    else Dy=0;
                    end
                    
                    D=[];Skdens=[];Laydist=[];
                    for bi=1:length(Bnd)-1
                       
                             if bi==length(Bnd)-1
                             indBnd=find(xi>=Bnd(bi));
                        else
                            indBnd=find(xi>=Bnd(bi)&xi<Bnd(bi+1));
                        end
                        
                       
                        
                        y1=nansum(avg(:,indBnd),2);
                        s1=skewness(y1);
                        Laydist=[Laydist y1];
                        
                        [mx ii]=max(y1);
                        jj=find(y1>mx*Perc);
                        
                        if ~isempty(jj)
                            D1=max(dx,abs(yi(jj(1))-yi(jj(end))));
                        else D1=0;
                        end
                        D=[D D1];
                        Skdens=[Skdens s1];
                    end
                    DW=[Dx ctextratio  Dy Dy./abs(ylim(1)-ylim(2)) D];
                    DistWidthCharge=[DistWidthCharge;DW];
                    Cellrecord{fn}.LaydistChgTp4Hp50=Laydist;                 
                    Skdenslayer=[Skdenslayer;Skdens];
                    Cellrecord{fn}.DistWidthdensTp4Hp50=DW;
                    Cellrecord{fn}.SkdenslayerTp4Hp50=Skdens;
                    
                end
            end
            
            
            if handles.avgTp==4&&abs(handles.avgHoldingPotential-10)<1&&isfield(Cellrecord{fn},'meanlatencyttxInd10')&&~isempty(Cellrecord{fn}.meanlatencyttxInd10)
                
                Soma= Cellrecord{fn}.SomaCoordinates;
                flipimg=Cellrecord{fn}.flipimg;
                pth=char(Cellrecord{fn}.Pth);
                stimcoordinates=Cellrecord{fn}.StimCoordinates;
                rotate_angle=Cellrecord{fn}.SpatialRotation;
                [stimCoordinates]=Matrixrotate_singlemap_GUI(Soma,stimcoordinates,flipimg,rotate_angle);
                X=stimCoordinates(1,:);
                Y=stimCoordinates(2,:);
                Zpk=zeros(size(X)).*NaN;
                Zchg=zeros(size(X)).*NaN;
                Zlat=zeros(size(X)).*NaN;
                if abs(Bnd(indb)-Bnd(1))>=handles.mindist&&abs(Bnd(indb)-Bnd(1))<=handles.maxdist
                     N=N+1;
                xaxis=[xaxis,X];
                yaxis=[yaxis,Y];
               
                dX=[dX,Cellrecord{fn}.PatternSpacing(1)];
                Zflag=Cellrecord{fn}.meanflagttxInd10;
                IndZ=find(Zflag>0);
                
                Zpk(IndZ)=Cellrecord{fn}.meanpeakttxInd10(IndZ);
                
                Zchg(IndZ)=Cellrecord{fn}.meanareattxInd10(IndZ);
                
                Zlat(IndZ)=Cellrecord{fn}.meanlatencyttxInd10(IndZ);
                Pchg=[];
                    Ppeak=[];
                    chg=[];
                    peak=[];
                    A10=[];
                     Mchg=[];
                    Mpeak=[];
                    for bi=1:length(Bnd)-1
                       
                            indBnd=find(X>=Bnd(bi)&X<Bnd(bi+1));
                     
                        
                        Pchg=[Pchg nansum(Zchg(indBnd))/nansum(Zchg)];
                        Ppeak=[Ppeak nansum(Zpk(indBnd))/nansum(Zpk)];
                        chg=[chg nansum(Zchg(indBnd))];
                        peak=[peak nansum(Zpk(indBnd))];
                         Mchg=[Mchg nanmean(Zchg(indBnd))];
                        Mpeak=[Mpeak nanmean(Zpk(indBnd))];
                        A10=[A10 nansum(Zpk(indBnd)>0).*Cellrecord{fn}.PatternSpacing(1).^2];
                    end
                Pchglayer=[Pchglayer;Pchg];
                Ppeaklayer=[Ppeaklayer;Ppeak];
                Pk=[Pk,Zpk];
               Chg=[Chg, Zchg];
                
                Lat=[Lat,Zlat];
                Den=[Den Zflag'];
                Cellrecord{fn}.PchgTp4Hp10=Pchg;
                      Cellrecord{fn}.PpeakTp4Hp10=Ppeak;
                      Cellrecord{fn}.chgTp4Hp10=chg;
                      Cellrecord{fn}.peakTp4Hp10=peak;
                       Cellrecord{fn}.AreaTp4Hp10=A10;
                       Cellrecord{fn}.MchgTp4Hp10=Mchg;
                      Cellrecord{fn}.MpeakTp4Hp10=Mpeak;
                    % Peak------------------------------------------
                    dx =Cellrecord{fn}.PatternSpacing(1); dy = dx;
                    xlim = [ceil(min(X))-dx-0.1 ceil(max(X))]; ylim = [ceil(min(Y))-dy-0.1 ceil(max(Y))];
                    
                    
                    xtick = xlim(1):dx:xlim(2); ytick = ylim(1):dy:ylim(2);
                    xi = xtick(1:end-1)+dx/2; yi = ytick(1:end-1)+dy/2;
                    avg = gridavg2(X,Y,Zpk,xtick,ytick,false,1,0);
                                       yy=nansum(avg,2);
                                       xx=nansum(avg);
                    [mxx ii]=max(xx);
                    jj=find(xx>mxx*Perc);
                    
                    if ~isempty(jj)
                        %                      ctext=abs(xi(jj(1))-xi(jj(end)));
                        ctextratio=abs(xi(jj(1))-xi(jj(end)))/abs(xlim(2)-xlim(1));
                        Dx=max(dx,abs(xi(jj(1))-xi(jj(end))));
                    else Dx=0;ctext=0;ctextratio=0;
                    end
                    [myy ii]=max(yy);
                    jj=find(yy>myy*Perc);
                    if ~isempty(jj)
                        Dy=max(dx,abs(yi(jj(1))-yi(jj(end))));
                    else Dy=0;
                    end
                  
                    D=[];Skpeak=[];Laydist=[];
                    for bi=1:length(Bnd)-1
                      if bi==length(Bnd)-1
                             indBnd=find(xi>=Bnd(bi));
                        else
                            indBnd=find(xi>=Bnd(bi)&xi<Bnd(bi+1));
                        end
                        
                    
                        
                        y1=nansum(avg(:,indBnd),2);
                        s1=skewness(y1);
                        Laydist=[Laydist y1];
                     
                        [mx ii]=max(y1);
                        jj=find(y1>mx*Perc);
                        
                        if ~isempty(jj)
                            D1=max(dx,abs(yi(jj(1))-yi(jj(end))));
                        else D1=0;
                        end
                        D=[D D1];
                        Skpeak=[Skpeak s1];
                    end
                   
                    DW=[Dx ctextratio  Dy Dy./abs(ylim(1)-ylim(2)) D];
                    DistWidthPeak=[DistWidthPeak;DW];
                    SkPeaklayer=[SkPeaklayer;Skpeak];
                    Cellrecord{fn}.LaydistPkTp4Hp10=Laydist;
                    Cellrecord{fn}.DistWidthPeakTp4Hp10=DW;
                    Cellrecord{fn}.SkPeaklayerTp4Hp10=Skpeak;
                    
                    
                    
                     %charge----------------------------------------------------------------------------------
                    avg = gridavg2(X,Y,Zchg,xtick,ytick,false,1,0);
                    xx=nansum(avg);
                    yy=nansum(avg,2);
                    [mxx ii]=max(xx);
                    jj=find(xx>mxx*Perc);
                    
                    if ~isempty(jj)
                        %                      ctext=abs(xi(jj(1))-xi(jj(end)));
                        ctextratio=abs(xi(jj(1))-xi(jj(end)))/abs(xlim(2)-xlim(1));
                        Dx=max(dx,abs(xi(jj(1))-xi(jj(end))));
                    else Dx=0;ctext=0;ctextratio=0;
                    end
                    
                    [myy ii]=max(yy);
                    jj=find(yy>myy*Perc);
                    if ~isempty(jj)
                        Dy=max(dx,abs(yi(jj(1))-yi(jj(end))));
                    else Dy=0;
                    end
             
                    D=[];Skchg=[];Laydist=[];
                    for bi=1:length(Bnd)-1
                      
                         if bi==length(Bnd)-1
                             indBnd=find(xi>=Bnd(bi));
                        else
                            indBnd=find(xi>=Bnd(bi)&xi<Bnd(bi+1));
                        end
                        
                       
                        y1=nansum(avg(:,indBnd),2);
                        s1=skewness(y1);
                        Laydist=[Laydist y1];
                   
                        [mx ii]=max(y1);
                        jj=find(y1>mx*Perc);
                        
                        if ~isempty(jj)
                            D1=max(dx,abs(yi(jj(1))-yi(jj(end))));
                        else D1=0;
                        end
                        D=[D D1];
                        Skchg=[Skchg s1];
                    end
                    
                    DW=[Dx ctextratio  Dy Dy./abs(ylim(1)-ylim(2)) D];
                    DistWidthCharge=[DistWidthCharge;DW];
                  
                    Skchglayer=[Skchglayer;Skchg];
                    Cellrecord{fn}.LaydistChgTp4Hp10=Laydist;                    
                    Cellrecord{fn}.DistWidthChargeTp4Hp10=DW;
                    Cellrecord{fn}.SkchglayerTp4Hp10=Skchg;
                    
                    
                    %density---------------------------------------------------------------
                    avg = gridavg2(X,Y,Zflag,xtick,ytick,false,1,0);
                    xx=nansum(avg);
                    yy=nansum(avg,2);
                    [mxx ii]=max(xx);
                    jj=find(xx>mxx*Perc);
                    
                    if ~isempty(jj)
                        %                      ctext=abs(xi(jj(1))-xi(jj(end)));
                        ctextratio=abs(xi(jj(1))-xi(jj(end)))/abs(xlim(2)-xlim(1));
                        Dx=max(dx,abs(xi(jj(1))-xi(jj(end))));
                    else Dx=0;ctext=0;ctextratio=0;
                    end
                    [myy ii]=max(yy);
                    jj=find(yy>myy*Perc);
                    if ~isempty(jj)
                        Dy=max(dx,abs(yi(jj(1))-yi(jj(end))));
                    else Dy=0;
                    end
                    
                    D=[];Skdens=[];Laydist=[];
                    for bi=1:length(Bnd)-1
                        
                     if bi==length(Bnd)-1
                             indBnd=find(xi>=Bnd(bi));
                        else
                            indBnd=find(xi>=Bnd(bi)&xi<Bnd(bi+1));
                        end
                        
                      
                        
                        y1=nansum(avg(:,indBnd),2);
                        s1=skewness(y1);
                        Laydist=[Laydist y1];
                        
                        [mx ii]=max(y1);
                        jj=find(y1>mx*Perc);
                        
                        if ~isempty(jj)
                            D1=max(dx,abs(yi(jj(1))-yi(jj(end))));
                        else D1=0;
                        end
                        D=[D D1];
                        Skdens=[Skdens s1];
                    end
                    DW=[Dx ctextratio  Dy Dy./abs(ylim(1)-ylim(2)) D];
                    DistWidthCharge=[DistWidthCharge;DW];
                    Cellrecord{fn}.LaydistChgTp4Hp10=Laydist;                 
                    Skdenslayer=[Skdenslayer;Skdens];
                    Cellrecord{fn}.DistWidthdensTp4Hp10=DW;
                    Cellrecord{fn}.SkdenslayerTp4Hp10=Skdens;
                    
                end
            end
            
            
            if handles.avgTp==5&&abs(handles.avgHoldingPotential+60)<1&&isfield(Cellrecord{fn},'meanlatencyotherInd60')&&~isempty(Cellrecord{fn}.meanlatencyotherInd60)
                
                Soma= Cellrecord{fn}.SomaCoordinates;
                flipimg=Cellrecord{fn}.flipimg;
                pth=char(Cellrecord{fn}.Pth)
                stimcoordinates=Cellrecord{fn}.StimCoordinates;
                rotate_angle=Cellrecord{fn}.SpatialRotation;
                [stimCoordinates]=Matrixrotate_singlemap_GUI(Soma,stimcoordinates,flipimg,rotate_angle);
                X=stimCoordinates(1,:);
                Y=stimCoordinates(2,:);
                
                Zpk=zeros(size(X)).*NaN;
                Zchg=zeros(size(X)).*NaN;
                Zlat=zeros(size(X)).*NaN;
                Ymat=zeros(size(X)).*NaN;
                Xmat=zeros(size(X)).*NaN;
                if abs(Bnd(indb)-Bnd(1))>=handles.mindist&&abs(Bnd(indb)-Bnd(1))<=handles.maxdist
                    N=N+1
      
                    xaxis=[xaxis,X];
                    yaxis=[yaxis,Y];
                    dX=[dX,Cellrecord{fn}.PatternSpacing(1)];  
                    Zflag=Cellrecord{fn}.meanflagotherInd60;
                    
                    IndZ=find(Zflag>0); 
                    dist=sqrt(X(IndZ).^2+Y(IndZ).^2);
                    dist=sort(dist);
                    if length(dist)
                    Dist8060=dist(ceil(length(dist)*0.8));
                    else Dist8060=0;
                    end
                    Zpk(IndZ)=Cellrecord{fn}.meanpeakotherInd60(IndZ);          
                    Zchg(IndZ)=Cellrecord{fn}.meanareaotherInd60(IndZ);                 
                    Zlat(IndZ)=Cellrecord{fn}.meanlatencyotherInd60(IndZ);
                    Ymat(IndZ)=Y(IndZ);
                         Xmat(IndZ)=X(IndZ);
                         
                          Xmat=Xmat./sqrt(Xmat.^2+Ymat.^2);
                                Ymat=Ymat./sqrt(Xmat.^2+Ymat.^2);
                    ZflagD=Cellrecord{fn}.meanflagDotherInd60;
                    IndZD=find(ZflagD>0); 
                    Cellrecord{fn}.DirectAreaTp5Hp60=length(IndZD).*Cellrecord{fn}.PatternSpacing(1).^2;
                    DirArea=[DirArea;length(IndZD).*Cellrecord{fn}.PatternSpacing(1).^2];
                    
                    Pchg=[];
                    Ppeak=[];
                    Parea=[];
                    chg=[];
                    peak=[];
                    Mchg=[];
                    Mpeak=[];
                     Medianchg=[];
                    Medianpeak=[];
                    A60=[];
                    AxisSum60=[];
                    expSum60=[];
                    ChgexpSum60=[];
                    PkexpSum60=[];
                    Dist8060Layer=[];
                    for bi=1:length(Bnd)-1
                       
                            indBnd=find(X>=Bnd(bi)&X<Bnd(bi+1));
                        
                        
                        Pchg=[Pchg nansum(Zchg(indBnd))/nansum(Zchg)];
                        Ppeak=[Ppeak nansum(Zpk(indBnd))/nansum(Zpk)];
                        Parea=[Parea nansum(Zflag(indBnd))/nansum(Zflag)];
                        chg=[chg nansum(Zchg(indBnd))];
                        peak=[peak nansum(Zpk(indBnd))];
                        Mchg=[Mchg nanmean(Zchg(indBnd))];
                        Mpeak=[Mpeak nanmean(Zpk(indBnd))];
                        Medianchg=[Medianchg nanmedian(Zchg(indBnd))];
                        Medianpeak=[Medianpeak nanmedian(Zpk(indBnd))];
                        A60=[A60 nansum(Zpk(indBnd)>0).*Cellrecord{fn}.PatternSpacing(1).^2];
                        AxisSum60=[AxisSum60 nansum(Ymat(indBnd))];
                        expSum60=[expSum60 abs(nansum((Xmat(indBnd)+i.*Ymat(indBnd))./sqrt(Xmat(indBnd).^2+Ymat(indBnd).^2)))./nansum(abs((Xmat(indBnd)+i.*Ymat(indBnd))./sqrt(Xmat(indBnd).^2+Ymat(indBnd).^2)))];
                        ChgexpSum60=[ChgexpSum60 abs(nansum(Zchg(indBnd).*(Xmat(indBnd)+i.*Ymat(indBnd))./sqrt(Xmat(indBnd).^2+Ymat(indBnd).^2)))./nansum(abs(Zchg(indBnd).*(Xmat(indBnd)+i.*Ymat(indBnd))./sqrt(Xmat(indBnd).^2+Ymat(indBnd).^2)))];
                        PkexpSum60=[PkexpSum60 abs(nansum(Zpk(indBnd).*(Xmat(indBnd)+i.*Ymat(indBnd))./sqrt(Xmat(indBnd).^2+Ymat(indBnd).^2)))./nansum(abs(Zpk(indBnd).*(Xmat(indBnd)+i.*Ymat(indBnd))./sqrt(Xmat(indBnd).^2+Ymat(indBnd).^2)))];
                        
                        indB=find(X(IndZ)>=Bnd(bi)&X(IndZ)<Bnd(bi+1));
                        if length(indB)>0
%                         dist=sqrt(X(IndZ(indB)).^2+Y(IndZ(indB)).^2);
                         dist=abs(Y(IndZ(indB)));
                        dist=sort(dist);
                         %Dist8060Layer=[Dist8060Layer abs(dist(ceil(length(dist)*0.1))-dist(ceil(length(dist)*0.9)))];
                        Dist8060Layer=[Dist8060Layer dist(ceil(length(dist)*0.8))];
                        else
                            Dist8060Layer=[Dist8060Layer 0];
                        end
                    end
                    
                    Pchglayer=[Pchglayer;Pchg];
                    Ppeaklayer=[Ppeaklayer;Ppeak];
                    Parealayer=[Parealayer;Parea];
                    
                    Pk=[Pk,Zpk];
                    Chg=[Chg, Zchg];
                    Lat=[Lat,Zlat];
                    Den=[Den Zflag'];
                    MCHG=[MCHG Mchg];
                    MPK=[MPK Mpeak];
                    
                    Cellrecord{fn}.PchgTp5Hp60=Pchg;
                    Cellrecord{fn}.PpeakTp5Hp60=Ppeak;
                    Cellrecord{fn}.PareaTp5Hp60=Parea;
                    Cellrecord{fn}.chgTp5Hp60=chg;
                    Cellrecord{fn}.peakTp5Hp60=peak;
                    Cellrecord{fn}.AreaTp5Hp60=A60;
                     Cellrecord{fn}.MchgTp5Hp60=Mchg;
                    Cellrecord{fn}.MpeakTp5Hp60=Mpeak;
                    Cellrecord{fn}.MedianchgTp5Hp60=Medianchg;
                    Cellrecord{fn}.MedianpeakTp5Hp60=Medianpeak;
                    Cellrecord{fn}.AxisSumTp5Hp60=AxisSum60;
                    Cellrecord{fn}.expSum60=expSum60;
                    Cellrecord{fn}.ChgexpSum60=ChgexpSum60;
                    Cellrecord{fn}.PkexpSum60=PkexpSum60;
                    Cellrecord{fn}.Dist8060Tp5Hp60=Dist8060;
                     Cellrecord{fn}.Dist8060LayerTp5Hp60=Dist8060Layer;
                    
                    % Peak---------------------------------------------------------------------------------------------------
                    dx =Cellrecord{fn}.PatternSpacing(1); dy = dx;
                    xlim = [ceil(min(X))-dx-0.1 ceil(max(X))]; ylim = [ceil(min(Y))-dy-0.1 ceil(max(Y))];
                    
                    
                    xtick = xlim(1):dx:xlim(2); ytick = ylim(1):dy:ylim(2);
                    xi = xtick(1:end-1)+dx/2; yi = ytick(1:end-1)+dy/2;
                    avg = gridavg2(X,Y,Zpk,xtick,ytick,false,1,0);
                    xx=nansum(avg);
                    yy=nansum(avg,2);
                    [mxx ii]=max(xx);
                    jj=find(xx>mxx*Perc);
                    
                    if ~isempty(jj)
                        %                      ctext=abs(xi(jj(1))-xi(jj(end)));
                        ctextratio=abs(xi(jj(1))-xi(jj(end)))/abs(max(X)-min(X));
                        Dx=max(dx,abs(xi(jj(1))-xi(jj(end))));
                    else Dx=0;ctext=0;ctextratio=0;
                    end
                    [myy ii]=max(yy);
                    jj=find(yy>myy*Perc);
                    if ~isempty(jj)
                        Dy=max(dx,abs(yi(jj(1))-yi(jj(end))));
                    else Dy=0;
                    end
                  
                    D=[];Skpeak=[];Laydist=[];
                    for bi=1:length(Bnd)-1
                        if bi==length(Bnd)-1
                             indBnd=find(xi>=Bnd(bi));
                        else
                            indBnd=find(xi>=Bnd(bi)&xi<Bnd(bi+1));
                        end
                        
                        y1=nansum(avg(:,indBnd),2);
                        s1=skewness(y1);
                        Laydist=[Laydist y1];
                     
                        [mx ii]=max(y1);
                        jj=find(y1>mx*Perc);
                        
                        if ~isempty(jj)
                            D1=max(dx,abs(yi(jj(1))-yi(jj(end))));
                        else D1=0;
                        end
                        D=[D D1];
                        Skpeak=[Skpeak s1];
                    end
                   
                    DW=[Dx ctextratio  Dy Dy./1000 D];
                    DistWidthPeak=[DistWidthPeak;DW];
                    SkPeaklayer=[SkPeaklayer;Skpeak];
                    Cellrecord{fn}.LaydistPkTp5Hp60=Laydist;
                    Cellrecord{fn}.DistWidthPeakTp5Hp60=DW;
                    Cellrecord{fn}.SkPeaklayerTp5Hp60=Skpeak;
                    
                    
                    %charge----------------------------------------------------------------------------------
                    avg = gridavg2(X,Y,Zchg,xtick,ytick,false,1,0);
                    xx=nansum(avg);
                    yy=nansum(avg,2);
                    [mxx ii]=max(xx);
                    jj=find(xx>mxx*Perc);
                    
                    if ~isempty(jj)
                        %                      ctext=abs(xi(jj(1))-xi(jj(end)));
                        ctextratio=abs(xi(jj(1))-xi(jj(end)))/abs(xlim(2)-xlim(1));
                        Dx=max(dx,abs(xi(jj(1))-xi(jj(end))));
                    else Dx=0;ctext=0;ctextratio=0;
                    end
                    
                    [myy ii]=max(yy);
                    jj=find(yy>myy*Perc);
                    if ~isempty(jj)
                        Dy=max(dx,abs(yi(jj(1))-yi(jj(end))));
                    else Dy=0;
                    end
             
                    D=[];Skchg=[];Laydist=[];
                    for bi=1:length(Bnd)-1
                        
                             if bi==length(Bnd)-1
                             indBnd=find(xi>=Bnd(bi));
                        else
                            indBnd=find(xi>=Bnd(bi)&xi<Bnd(bi+1));
                        end
                        
                        
                        y1=nansum(avg(:,indBnd),2);
                        s1=skewness(y1);
                        Laydist=[Laydist y1];
                   
                        [mx ii]=max(y1);
                        jj=find(y1>mx*Perc);
                        
                        if ~isempty(jj)
                            D1=max(dx,abs(yi(jj(1))-yi(jj(end))));
                        else D1=0;
                        end
                        D=[D D1];
                        Skchg=[Skchg s1];
                    end
                    
                    DW=[Dx ctextratio  Dy Dy./abs(1000) D];
                    DistWidthCharge=[DistWidthCharge;DW];
                  
                    Skchglayer=[Skchglayer;Skchg];
                    Cellrecord{fn}.LaydistChgTp5Hp60=Laydist;
                    
                    Cellrecord{fn}.DistWidthChargeTp5Hp60=DW;
                    Cellrecord{fn}.SkchglayerTp5Hp60=Skchg;
                    
                    
                    %density---------------------------------------------------------------
                    avg = gridavg2(X,Y,Zflag,xtick,ytick,false,1,0);
                    xx=nansum(avg);
                    yy=nansum(avg,2);
                    [mxx ii]=max(xx);
                    jj=find(xx>mxx*Perc);
                    
                    if ~isempty(jj)
                        %                      ctext=abs(xi(jj(1))-xi(jj(end)));
                        ctextratio=abs(xi(jj(1))-xi(jj(end)))/abs(xlim(2)-xlim(1));
                        Dx=max(dx,abs(xi(jj(1))-xi(jj(end))));
                    else Dx=0;ctext=0;ctextratio=0;
                    end
                    [myy ii]=max(yy);
                    jj=find(yy>myy*Perc);
                    if ~isempty(jj)
                        Dy=max(dx,abs(yi(jj(1))-yi(jj(end))));
                    else Dy=0;
                    end
                    
                    D=[];Skdens=[];Laydist=[];
                    for bi=1:length(Bnd)-1
                       
                             if bi==length(Bnd)-1
                             indBnd=find(xi>=Bnd(bi));
                        else
                            indBnd=find(xi>=Bnd(bi)&xi<Bnd(bi+1));
                        end
                        
                        
                        y1=nansum(avg(:,indBnd),2);
                        s1=skewness(y1);
                        Laydist=[Laydist y1];
                        
                        [mx ii]=max(y1);
                        jj=find(y1>mx*Perc);
                        
                        if ~isempty(jj)
                            D1=max(dx,abs(yi(jj(1))-yi(jj(end))));
                        else D1=0;
                        end
                        D=[D D1];
                        Skdens=[Skdens s1];
                    end
                    DW=[Dx ctextratio  Dy Dy./abs(1000) D];
                    DistWidthDense=[DistWidthDense;DW];
                    Cellrecord{fn}.LaydistDensTp5Hp60=Laydist;                 
                    Skdenslayer=[Skdenslayer;Skdens];
                    Cellrecord{fn}.DistWidthdensTp5Hp60=DW;
                    Cellrecord{fn}.SkdenslayerTp5Hp60=Skdens;
                    
                end
            end
            
            if handles.avgTp==5&&abs(handles.avgHoldingPotential-10)<1&&isfield(Cellrecord{fn},'meanlatencyotherInd10')&&~isempty(Cellrecord{fn}.meanlatencyotherInd10)
                
                Soma= Cellrecord{fn}.SomaCoordinates;
                flipimg=Cellrecord{fn}.flipimg;
                pth=char(Cellrecord{fn}.Pth);
                stimcoordinates=Cellrecord{fn}.StimCoordinates;
                rotate_angle=Cellrecord{fn}.SpatialRotation;
                [stimCoordinates]=Matrixrotate_singlemap_GUI(Soma,stimcoordinates,flipimg,rotate_angle);
                X=stimCoordinates(1,:);
                Y=stimCoordinates(2,:);
                Zpk=zeros(size(X)).*NaN;
                Zchg=zeros(size(X)).*NaN;
                Zlat=zeros(size(X)).*NaN;
                  Ymat=zeros(size(X)).*NaN;
                   Xmat=zeros(size(X)).*NaN;
                if abs(Bnd(indb)-Bnd(1))>=handles.mindist&&abs(Bnd(indb)-Bnd(1))<=handles.maxdist
                    N=N+1;
                    xaxis=[xaxis,X];
                    yaxis=[yaxis,Y];
                    
                    dX=[dX,Cellrecord{fn}.PatternSpacing(1)];
                    Zflag=Cellrecord{fn}.meanflagotherInd10;
                    IndZ=find(Zflag>0);
                    if length(IndZ)
                     dist=sqrt(X(IndZ).^2+Y(IndZ).^2);
                    dist=sort(dist);
                    Dist8010=dist(ceil(length(dist)*0.8));
                    else Dist8010=0;
                    end
                    Zpk(IndZ)=Cellrecord{fn}.meanpeakotherInd10(IndZ);
                    
                    Zchg(IndZ)=Cellrecord{fn}.meanareaotherInd10(IndZ);
                    
                    Zlat(IndZ)=Cellrecord{fn}.meanlatencyotherInd10(IndZ);
                    Ymat(IndZ)=Y(IndZ);
                    Xmat(IndZ)=X(IndZ);
                    Pchg=[];
                    Ppeak=[];
                    Parea=[];
                    chg=[];
                    peak=[];
                    A10=[];
                       Mchg=[];
                    Mpeak=[];
                          Medianchg=[];
                    Medianpeak=[];
                     AxisSum10=[];
                     expSum10=[];
                    ChgexpSum10=[];
                    PkexpSum10=[];
                    
                    Dist8010Layer=[];
                    Xmat=Xmat./sqrt(Xmat.^2+Ymat.^2);
                                Ymat=Ymat./sqrt(Xmat.^2+Ymat.^2);
                    for bi=1:length(Bnd)-1
                        
                            indBnd=find(X>=Bnd(bi)&X<Bnd(bi+1));
                      
                        
                        Pchg=[Pchg nansum(Zchg(indBnd))/nansum(Zchg)];
                        Ppeak=[Ppeak nansum(Zpk(indBnd))/nansum(Zpk)];
                        Parea=[Parea nansum(Zflag(indBnd))/nansum(Zflag)];
                       chg=[chg nansum(Zchg(indBnd))];
                        peak=[peak nansum(Zpk(indBnd))];
                        Mchg=[Mchg nanmean(Zchg(indBnd))];
                        Mpeak=[Mpeak nanmean(Zpk(indBnd))];
                           Medianchg=[Medianchg nanmedian(Zchg(indBnd))];
                        Medianpeak=[Medianpeak nanmedian(Zpk(indBnd))];
                        A10=[A10 nansum(Zpk(indBnd)>0).*Cellrecord{fn}.PatternSpacing(1).^2];
                         AxisSum10=[AxisSum10 nansum(Ymat(indBnd))];
                         expSum10=[expSum10 abs(nansum((Xmat(indBnd)+i.*Ymat(indBnd))./sqrt(Xmat(indBnd).^2+Ymat(indBnd).^2)))./nansum(abs((Xmat(indBnd)+i.*Ymat(indBnd))./sqrt(Xmat(indBnd).^2+Ymat(indBnd).^2)))];
                        ChgexpSum10=[ChgexpSum10 abs(nansum(Zchg(indBnd).*(Xmat(indBnd)+i.*Ymat(indBnd))./sqrt(Xmat(indBnd).^2+Ymat(indBnd).^2)))./nansum(abs(Zchg(indBnd).*(Xmat(indBnd)+i.*Ymat(indBnd))./sqrt(Xmat(indBnd).^2+Ymat(indBnd).^2)))];
                        PkexpSum10=[PkexpSum10 abs(nansum(Zpk(indBnd).*(Xmat(indBnd)+i.*Ymat(indBnd))./sqrt(Xmat(indBnd).^2+Ymat(indBnd).^2)))./nansum(abs(Zpk(indBnd).*(Xmat(indBnd)+i.*Ymat(indBnd))./sqrt(Xmat(indBnd).^2+Ymat(indBnd).^2)))];
                         indB=find(X(IndZ)>=Bnd(bi)&X(IndZ)<Bnd(bi+1));
                         if length(indB)>0
%                         dist=sqrt(X(IndZ(indB)).^2+Y(IndZ(indB)).^2);
                   dist=abs(Y(IndZ(indB)));
                        dist=sort(dist);
                         %Dist8010Layer=[Dist8010Layer abs(dist(ceil(length(dist)*0.1))-dist(ceil(length(dist)*0.9)))];
                          Dist8010Layer=[Dist8010Layer dist(ceil(length(dist)*0.8))];
                         else
                             Dist8010Layer=[Dist8010Layer 0];
                         end
                    end
                    Pchglayer=[Pchglayer;Pchg];
                    Ppeaklayer=[Ppeaklayer;Ppeak];
                    Parealayer=[Parealayer;Parea];
                    Pk=[Pk,Zpk];
                    Chg=[Chg, Zchg];                    
                    Lat=[Lat,Zlat];
                    Den=[Den Zflag'];
                    Cellrecord{fn}.PchgTp5Hp10=Pchg;
                    Cellrecord{fn}.PpeakTp5Hp10=Ppeak;
                    Cellrecord{fn}.PareaTp5Hp10=Parea;
                    Cellrecord{fn}.chgTp5Hp10=chg;
                    Cellrecord{fn}.peakTp5Hp10=peak;
                    Cellrecord{fn}.AreaTp5Hp10=A10;
                    Cellrecord{fn}.MchgTp5Hp10=Mchg;
                    Cellrecord{fn}.MpeakTp5Hp10=Mpeak;
                    Cellrecord{fn}.AxisSumTp5Hp10=AxisSum10;
                     
                    Cellrecord{fn}.expSum10=expSum10;
                    Cellrecord{fn}.ChgexpSum10=ChgexpSum10;
                    Cellrecord{fn}.PkexpSum10=PkexpSum10;
                             Cellrecord{fn}.MedianchgTp5Hp10=Medianchg;
                    Cellrecord{fn}.MedianpeakTp5Hp10=Medianpeak;
                     Cellrecord{fn}.Dist8010Tp5Hp60=Dist8010;
                     Cellrecord{fn}.Dist8010LayerTp5Hp10=Dist8010Layer;
                    
                    % Peak------------------------------------------
                    dx =Cellrecord{fn}.PatternSpacing(1); dy = dx;
                    xlim = [ceil(min(X))-dx-0.1 ceil(max(X))]; ylim = [ceil(min(Y))-dy-0.1 ceil(max(Y))];
                    
                    
                    xtick = xlim(1):dx:xlim(2); ytick = ylim(1):dy:ylim(2);
                    xi = xtick(1:end-1)+dx/2; yi = ytick(1:end-1)+dy/2;
                    avg = gridavg2(X,Y,Zpk,xtick,ytick,false,1,0);
                                       yy=nansum(avg,2);
                                       xx=nansum(avg);
                    [mxx ii]=max(xx);
                    jj=find(xx>mxx*Perc);
                    
                    if ~isempty(jj)
                        %                      ctext=abs(xi(jj(1))-xi(jj(end)));
                        ctextratio=abs(xi(jj(1))-xi(jj(end)))/abs(xlim(2)-xlim(1));
                        Dx=max(dx,abs(xi(jj(1))-xi(jj(end))));
                    else Dx=0;ctext=0;ctextratio=0;
                    end
                    [myy ii]=max(yy);
                    jj=find(yy>myy*Perc);
                    if ~isempty(jj)
                        Dy=max(dx,abs(yi(jj(1))-yi(jj(end))));
                    else Dy=0;
                    end
                  
                    D=[];Skpeak=[];Laydist=[];
                    for bi=1:length(Bnd)-1
                        
                          if bi==length(Bnd)-1
                             indBnd=find(xi>=Bnd(bi));
                        else
                            indBnd=find(xi>=Bnd(bi)&xi<Bnd(bi+1));
                        end
                        
                        
                        y1=nansum(avg(:,indBnd),2);
                        s1=skewness(y1);
                        Laydist=[Laydist y1];
                     
                        [mx ii]=max(y1);
                        jj=find(y1>mx*Perc);
                        
                        if ~isempty(jj)
                            D1=max(dx,abs(yi(jj(1))-yi(jj(end))));
                        else D1=0;
                        end
                        D=[D D1];
                        Skpeak=[Skpeak s1];
                    end
                   
                    DW=[Dx ctextratio  Dy Dy./abs(ylim(1)-ylim(2)) D];
                    DistWidthPeak=[DistWidthPeak;DW];
                    SkPeaklayer=[SkPeaklayer;Skpeak];
                    Cellrecord{fn}.LaydistPkTp5Hp10=Laydist;
                    Cellrecord{fn}.DistWidthPeakTp5Hp10=DW;
                    Cellrecord{fn}.SkPeaklayerTp5Hp10=Skpeak;
                    
                    
                    
                     %charge----------------------------------------------------------------------------------
                    avg = gridavg2(X,Y,Zchg,xtick,ytick,false,1,0);
                    xx=nansum(avg);
                    yy=nansum(avg,2);
                    [mxx ii]=max(xx);
                    jj=find(xx>mxx*Perc);
                    
                    if ~isempty(jj)
                        %                      ctext=abs(xi(jj(1))-xi(jj(end)));
                        ctextratio=abs(xi(jj(1))-xi(jj(end)))/abs(xlim(2)-xlim(1));
                        Dx=max(dx,abs(xi(jj(1))-xi(jj(end))));
                    else Dx=0;ctext=0;ctextratio=0;
                    end
                    
                    [myy ii]=max(yy);
                    jj=find(yy>myy*Perc);
                    if ~isempty(jj)
                        Dy=max(dx,abs(yi(jj(1))-yi(jj(end))));
                    else Dy=0;
                    end
             
                    D=[];Skchg=[];Laydist=[];
                    for bi=1:length(Bnd)-1
                       
                         if bi==length(Bnd)-1
                             indBnd=find(xi>=Bnd(bi));
                        else
                            indBnd=find(xi>=Bnd(bi)&xi<Bnd(bi+1));
                        end
                        
                        
                        
                        y1=nansum(avg(:,indBnd),2);
                        s1=skewness(y1);
                        Laydist=[Laydist y1];
                   
                        [mx ii]=max(y1);
                        jj=find(y1>mx*Perc);
                        
                        if ~isempty(jj)
                            D1=max(dx,abs(yi(jj(1))-yi(jj(end))));
                        else D1=0;
                        end
                        D=[D D1];
                        Skchg=[Skchg s1];
                    end
                    
                    DW=[Dx ctextratio  Dy Dy./abs(ylim(1)-ylim(2)) D];
                    DistWidthCharge=[DistWidthCharge;DW];
                  
                    Skchglayer=[Skchglayer;Skchg];
                    Cellrecord{fn}.LaydistChgTp5Hp10=Laydist;                    
                    Cellrecord{fn}.DistWidthChargeTp5Hp10=DW;
                    Cellrecord{fn}.SkchglayerTp5Hp10=Skchg;
                    
                    
                    %density---------------------------------------------------------------
                    avg = gridavg2(X,Y,Zflag,xtick,ytick,false,1,0);
                    xx=nansum(avg);
                    yy=nansum(avg,2);
                    [mxx ii]=max(xx);
                    jj=find(xx>mxx*Perc);
                    
                    if ~isempty(jj)
                        %                      ctext=abs(xi(jj(1))-xi(jj(end)));
                        ctextratio=abs(xi(jj(1))-xi(jj(end)))/abs(xlim(2)-xlim(1));
                        Dx=max(dx,abs(xi(jj(1))-xi(jj(end))));
                    else Dx=0;ctext=0;ctextratio=0;
                    end
                    [myy ii]=max(yy);
                    jj=find(yy>myy*Perc);
                    if ~isempty(jj)
                        Dy=max(dx,abs(yi(jj(1))-yi(jj(end))));
                    else Dy=0;
                    end
                    
                    D=[];Skdens=[];Laydist=[];
                    for bi=1:length(Bnd)-1
                        
                            if bi==length(Bnd)-1
                             indBnd=find(xi>=Bnd(bi));
                        else
                            indBnd=find(xi>=Bnd(bi)&xi<Bnd(bi+1));
                        end
                        
                     
                        
                        y1=nansum(avg(:,indBnd),2);
                        s1=skewness(y1);
                        Laydist=[Laydist y1];
                        
                        [mx ii]=max(y1);
                        jj=find(y1>mx*Perc);
                        
                        if ~isempty(jj)
                            D1=max(dx,abs(yi(jj(1))-yi(jj(end))));
                        else D1=0;
                        end
                        D=[D D1];
                        Skdens=[Skdens s1];
                    end
                    DW=[Dx ctextratio  Dy Dy./abs(ylim(1)-ylim(2)) D];
                    DistWidthDense=[DistWidthDense;DW];
                    Cellrecord{fn}.LaydistDensTp5Hp10=Laydist;                 
                    Skdenslayer=[Skdenslayer;Skdens];
                    Cellrecord{fn}.DistWidthdensTp5Hp10=DW;
                    Cellrecord{fn}.SkdenslayerTp5Hp10=Skdens;
                    
                end
            end
            
            
        end
      
      
    end
      Num_cell=Cg.Num_cell;
      filepath=Cg.filepath;
    save(fullfile(folderMeanpath,a(i).name),'Cellrecord','Num_cell','filepath'); 
    end
end
end
  Cellsum.avgPk=Pk;
   Cellsum.avgChg=Chg;
   Cellsum.Lat=Lat;
  Cellsum.dX=dX;
  Cellsum.AvgAge=AGE;
  Cellsum.Bnd=BD;
  Cellsum.avgxaxis=xaxis;
  Cellsum.avgyaxis=yaxis;
  Cellsum.avgDensity=Den;
  Cellsum.avgCellNum=N;
  Cellsum.Pchglayer=Pchglayer;
  Cellsum.Ppeaklayer=Ppeaklayer;
  Cellsum.DistWidthCharge=DistWidthCharge;
  Cellsum.Skchglayer=Skchglayer;
  Cellsum.DistWidthPeak=DistWidthPeak;
  Cellsum.SkPeaklayer=SkPeaklayer;
  Cellsum.DirArea=DirArea;
  %silent
   Cellsum.avgPksi=Pksi;
  Cellsum.avgChgsi=Chgsi;
  Cellsum.Latsi=Latsi;
  Cellsum.dXsi=dXsi;
  Cellsum.AvgAgesi=AGEsi;
  Cellsum.Bndsi=BDsi;
  Cellsum.avgxaxissi=xaxissi;
  Cellsum.avgyaxissi=yaxissi;
  Cellsum.avgDensitysi=Densi;
  Cellsum.avgCellNumsi=Nsi;
  Cellsum.Pchgsilayer=Pchgsilayer;
  Cellsum.Ppeaksilayer=Ppeaksilayer;
  Cellsum.DistWidthChargesi=DistWidthChargesi;
  Cellsum.DistWidthPeaksi=DistWidthPeaksi;
  Cellsum.ratio=AR;
  eval(['Cellsum' '.' sprintf('PchglayerTp%iHP%iDrexp%iAge%ito%i',handles.avgTp,abs(handles.avgHoldingPotential),handles.avgdarkexp,handles.minage,handles.maxage) '=Pchglayer' ';']);
   eval([sprintf('Cellsum.PpeaklayerTp%iHP%iDrexp%iAge%ito%i',handles.avgTp,abs(handles.avgHoldingPotential),handles.avgdarkexp,handles.minage,handles.maxage) '=Ppeaklayer' ';']);
   eval([sprintf('Cellsum.DistWidthChargeTp%iHP%iDrexp%iAge%ito%i',handles.avgTp,abs(handles.avgHoldingPotential),handles.avgdarkexp,handles.minage,handles.maxage) '=DistWidthCharge' ';']);
   eval([sprintf('Cellsum.DistWidthPeakTp%iHP%iDrexp%iAge%ito%i',handles.avgTp,abs(handles.avgHoldingPotential),handles.avgdarkexp,handles.minage,handles.maxage) '=DistWidthPeak' ';']);
   eval([sprintf('Cellsum.SkchglayerTp%iHP%iDrexp%iAge%ito%i',handles.avgTp,abs(handles.avgHoldingPotential),handles.avgdarkexp,handles.minage,handles.maxage) '=DistWidthCharge' ';']);
   eval([sprintf('Cellsum.SkPeaklayerTp%iHP%iDrexp%iAge%ito%i',handles.avgTp,abs(handles.avgHoldingPotential),handles.avgdarkexp,handles.minage,handles.maxage) '=DistWidthPeak' ';']);
   eval([sprintf('Cellsum.avgxaxisTp%iHP%iDrexp%iAge%ito%i',handles.avgTp,abs(handles.avgHoldingPotential),handles.avgdarkexp,handles.minage,handles.maxage) '=xaxis' ';']);
   eval([sprintf('Cellsum.avgyaxisTp%iHP%iDrexp%iAge%ito%i',handles.avgTp,abs(handles.avgHoldingPotential),handles.avgdarkexp,handles.minage,handles.maxage) '=yaxis' ';']);
   eval([sprintf('Cellsum.avgPkTp%iHP%iDrexp%iAge%ito%i',handles.avgTp,abs(handles.avgHoldingPotential),handles.avgdarkexp,handles.minage,handles.maxage) '=Pk' ';']);
   eval([sprintf('Cellsum.avgChgTp%iHP%iDrexp%iAge%ito%i',handles.avgTp,abs(handles.avgHoldingPotential),handles.avgdarkexp,handles.minage,handles.maxage) '=Chg' ';']);
   eval([sprintf('Cellsum.avgdensityTp%iHP%iDrexp%iAge%ito%i',handles.avgTp,abs(handles.avgHoldingPotential),handles.avgdarkexp,handles.minage,handles.maxage) '=Den' ';']);
   eval([sprintf('Cellsum.avgDirAreaTp%iHP%iDrexp%iAge%ito%i',handles.avgTp,abs(handles.avgHoldingPotential),handles.avgdarkexp,handles.minage,handles.maxage) '=DirArea' ';']);
   % for silent
   if handles.avgTp==2
    eval(['Cellsum' '.' sprintf('PchglayerTp%iHP%iDrexp%isiAge%ito%i',handles.avgTp,abs(handles.avgHoldingPotential),handles.avgdarkexp,handles.minage,handles.maxage) '=Pchgsilayer' ';']);
   eval([sprintf('Cellsum.PpeaklayerTp%iHP%iDrexp%isiAge%ito%i',handles.avgTp,abs(handles.avgHoldingPotential),handles.avgdarkexp,handles.minage,handles.maxage) '=Ppeaksilayer' ';']);
   eval([sprintf('Cellsum.DistWidthChargeTp%iHP%iDrexp%isiAge%ito%i',handles.avgTp,abs(handles.avgHoldingPotential),handles.avgdarkexp,handles.minage,handles.maxage) '=DistWidthChargesi' ';']);
   eval([sprintf('Cellsum.DistWidthPeakTp%iHP%iDrexp%isiAge%ito%i',handles.avgTp,abs(handles.avgHoldingPotential),handles.avgdarkexp,handles.minage,handles.maxage) '=DistWidthPeaksi' ';']);
   eval([sprintf('Cellsum.SkchglayerTp%iHP%iDrexp%isiAge%ito%i',handles.avgTp,abs(handles.avgHoldingPotential),handles.avgdarkexp,handles.minage,handles.maxage) '=DistWidthChargesi' ';']);
   eval([sprintf('Cellsum.SkPeaklayerTp%iHP%iDrexp%isiAge%ito%i',handles.avgTp,abs(handles.avgHoldingPotential),handles.avgdarkexp,handles.minage,handles.maxage) '=DistWidthPeaksi' ';']);
   eval([sprintf('Cellsum.avgxaxisTp%iHP%iDrexp%isiAge%ito%i',handles.avgTp,abs(handles.avgHoldingPotential),handles.avgdarkexp,handles.minage,handles.maxage) '=xaxissi' ';']);
   eval([sprintf('Cellsum.avgyaxisTp%iHP%iDrexp%isiAge%ito%i',handles.avgTp,abs(handles.avgHoldingPotential),handles.avgdarkexp,handles.minage,handles.maxage) '=yaxissi' ';']);
   eval([sprintf('Cellsum.avgPkTp%iHP%iDrexp%isiAge%ito%i',handles.avgTp,abs(handles.avgHoldingPotential),handles.avgdarkexp,handles.minage,handles.maxage) '=Pksi' ';']);
   eval([sprintf('Cellsum.avgChgTp%iHP%iDrexp%isiAge%ito%i',handles.avgTp,abs(handles.avgHoldingPotential),handles.avgdarkexp,handles.minage,handles.maxage) '=Chgsi' ';']);
   eval([sprintf('Cellsum.avgdensityTp%iHP%iDrexp%isiAge%ito%i',handles.avgTp,abs(handles.avgHoldingPotential),handles.avgdarkexp,handles.minage,handles.maxage) '=Densi' ';']);
   eval([sprintf('Cellsum.avgAreaRatioTp%iHP%iDrexp%isiAge%ito%i',handles.avgTp,abs(handles.avgHoldingPotential),handles.avgdarkexp,handles.minage,handles.maxage) '=AR' ';']);
   end
   filesum=sprintf('SilentsynapseSummaryAge%ito%i.mat',handles.minage,handles.maxage)
  save(fullfile(folderMeanpath,filesum),'Cellsum') 
  if N<1
     set(handles.error_run,'string','No maps meets your standards')
        else
        set(handles.error_run,'string','')
  end
end

