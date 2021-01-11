function [ output_args ] = COMPARE_Maps_3(Control_path, Object_path,savepath)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%type=1; %use ranksum
folderMeanpath{1}=Control_path;
folderMeanpath{2}=Object_path;

close all;



 N=[];
for cmpind=1:2
    a=dir(fullfile(folderMeanpath{cmpind},'*.mat'));
    Cellfile=size(a);
    
    DirArea=[];
    Ag=[];
    Bnd=[];
    n=0;
    Nc=(Cellfile(1,1))*5;
    
    Xaxis=cell(Nc,1);
    Yaxis=cell(Nc,1);
    Zpk60avg=cell(Nc,1);
    Zchg60avg=cell(Nc,1);
    Zden60avg=cell(Nc,1);
    Zpk10avg=cell(Nc,1);
    Zchg10avg=cell(Nc,1);
    Zden10avg=cell(Nc,1);
    Zlat10avg=cell(Nc,1);
    Zlat60avg=cell(Nc,1);
    
    ZchgD60avg=cell(Nc,1);
    ZdenD60avg=cell(Nc,1);
    
    ZlatD60avg=cell(Nc,1);
    for i= 1:1:Cellfile
        Cg= load(fullfile(folderMeanpath{cmpind},a(i).name));
        
        
        if isfield(Cg,'Cellrecord')
            Cellrecord=Cg.Cellrecord;
            for fn=1:Cg.Num_cell+1
                
                if isfield(Cellrecord{fn},'DistWidthPeakTp1Hp60')||isfield(Cellrecord{fn},'DistWidthPeakTp1Hp10')
                    if isfield(Cellrecord{fn},'DistWidthPeakTp1Hp60')
                        NL=length(Cellrecord{fn}.DistWidthPeakTp1Hp60);
                    else NL=length(Cellrecord{fn}.DistWidthPeakTp1Hp10);
                    end
                    Ag=[Ag;Cellrecord{fn}.age];
                    Bnd=[Bnd;Cellrecord{fn}.Boundry];
                    n=n+1;
                    
                    Soma= Cellrecord{fn}.SomaCoordinates;
                    flipimg=Cellrecord{fn}.flipimg;
                    pth=char(Cellrecord{fn}.Pth);
                    stimcoordinates=Cellrecord{fn}.StimCoordinates;
                    rotate_angle=Cellrecord{fn}.SpatialRotation;
                    [stimCoordinates]=Matrixrotate_singlemap_GUI(Soma,stimcoordinates,flipimg,rotate_angle);
                    X=stimCoordinates(1,:);
                    Y=stimCoordinates(2,:);
                    Xaxis{n}=X;
                    Yaxis{n}=Y;
                    Zpk60=zeros(size(X)).*NaN;
                    Zchg60=zeros(size(X)).*NaN;
                    Zden60=zeros(size(X));
                    Zpk10=zeros(size(X)).*NaN;
                    Zden10=zeros(size(X));
                    Zchg10=zeros(size(X)).*NaN;
                    Zlat60=zeros(size(X)).*NaN;
                    Zlat10=zeros(size(X)).*NaN;
                    ZlatD60=zeros(size(X)).*NaN;
                    ZchgD60=zeros(size(X)).*NaN;
                    ZdenD60=zeros(size(X));
                    
                    ZlatD10=zeros(size(X)).*NaN;
                    
                    if isfield(Cellrecord{fn},'DistWidthPeakTp1Hp60')
                        
                        LdistPk60{n}=Cellrecord{fn}.LaydistPkTp1Hp60;
                        LdistChg60{n}=Cellrecord{fn}.LaydistChgTp1Hp60;
                        
                        Indden=find((Cellrecord{fn}.meanflagDHighMg60>0));
                        DirArea=[DirArea,length(Indden).*Cellrecord{fn}.PatternSpacing(1).^2];
                        Indpk=find(~isnan(Cellrecord{fn}.meanpeakHighMg60));
                        Zpk60(Indpk)=Cellrecord{fn}.meanpeakHighMg60(Indpk);
                        Indchg=find(~isnan(Cellrecord{fn}.meanareaHighMg60));
                        Zchg60(Indchg)=Cellrecord{fn}.meanareaHighMg60(Indchg);
                        Indlat=find(~isnan(Cellrecord{fn}.meanlatencyHighMg60));
                        Zlat60(Indlat)=Cellrecord{fn}.meanlatencyHighMg60(Indlat);
                        Indden=find((Cellrecord{fn}.meanflagHighMg60>0));
                        Zden60(Indden)=Cellrecord{fn}.meanflagHighMg60(Indden);
                        
                        
                        Indchg=find(~isnan(Cellrecord{fn}.meanDirareaHighMg60));
                        ZchgD60(Indchg)=Cellrecord{fn}.meanDirareaHighMg60(Indchg);
                        Indlat=find(~isnan(Cellrecord{fn}.meanDirlatencyHighMg60));
                        ZlatD60(Indlat)=Cellrecord{fn}.meanDirlatencyHighMg60(Indlat);
                        Indden=find((Cellrecord{fn}.meanflagDHighMg60>0));
                        ZdenD60(Indden)=Cellrecord{fn}.meanflagDHighMg60(Indden);
                        
                        
                        
                    else Nanmaxtric= ones(1,NL).*NaN;
                        
                        LdistPk60{n}=NaN;
                        LdistChg60{n}=NaN;
                    end
                    if isfield(Cellrecord{fn},'DistWidthPeakTp1Hp10')
                        
                        LdistPk10{n}=Cellrecord{fn}.LaydistPkTp1Hp10;
                        LdistChg10{n}=Cellrecord{fn}.LaydistChgTp1Hp10;
                        Indpk=find(~isnan(Cellrecord{fn}.meanpeakHighMg10));
                        Zpk10(Indpk)=Cellrecord{fn}.meanpeakHighMg10(Indpk);
                        Indchg=find(~isnan(Cellrecord{fn}.meanareaHighMg10));
                        Zchg10(Indchg)=Cellrecord{fn}.meanareaHighMg10(Indchg);
                        Indlat=find(~isnan(Cellrecord{fn}.meanlatencyHighMg10));
                        Zlat10(Indlat)=Cellrecord{fn}.meanlatencyHighMg10(Indlat);
                        Indden=find((Cellrecord{fn}.meanflagHighMg10>0));
                        Zden10(Indden)=Cellrecord{fn}.meanflagHighMg10(Indden);
                        
                        
                    else
                        
                        LdistPk10{n}=NaN;
                        LdistChg10{n}=NaN;
                    end
                    
                    Zpk60avg{n}=Zpk60;
                    Zchg60avg{n}=Zchg60;
                    Zden60avg{n}=Zden60;
                    Zpk10avg{n}=Zpk10;
                    Zchg10avg{n}=Zchg10;
                    Zden10avg{n}=Zden10;
                    Zlat10avg{n}=Zlat10;
                    Zlat60avg{n}=Zlat60;
                    
                    ZchgD60avg{n}=ZchgD60;
                    ZdenD60avg{n}=ZdenD60;
                    
                    ZlatD60avg{n}=ZlatD60;
                    
                end
                
                
                
            end
        end
    end
    Cs{cmpind}.Zpk60avg=Zpk60avg;
    Cs{cmpind}.Zchg60avg=Zchg60avg;
    Cs{cmpind}.Zden60avg=Zden60avg;
    Cs{cmpind}.Zpk10avg=Zpk10avg;
    
    Cs{cmpind}.Zchg10avg=Zchg10avg;
    Cs{cmpind}.Zden10avg=Zden10avg;
    Cs{cmpind}.Zlat10avg=Zlat10avg;
    Cs{cmpind}. Zlat60avg= Zlat60avg;
    
    
    Cs{cmpind}.ZchgD60avg=ZchgD60avg;
    Cs{cmpind}.ZdenD60avg=ZdenD60avg;
    
    Cs{cmpind}.ZlatD60avg=ZlatD60avg;
    Cs{cmpind}.Yaxis=Yaxis;
    Cs{cmpind}.DirArea=DirArea;
    
    Cs{cmpind}.Xaxis=Xaxis;
    Cs{cmpind}.Bnd=Bnd;
    
    Cs{cmpind}.LdistPk60=LdistPk60;
    Cs{cmpind}.LdistPk10=LdistPk10;
    Cs{cmpind}.LdistChg60=LdistChg60;
    Cs{cmpind}.LdistChg10=LdistChg10;
    N=[N;n]
    
    
end
 dx=30;
    dy=30;
    
xlim = [-420-0.1 960];
    ylim = [-480-0.1 480+0.1];
    
    xtick = xlim(1):dx:xlim(2);
    ytick = ylim(1):dy:ylim(2);
    xi = xtick(1:end-1)+dx/2;
    yi = ytick(1:end-1)+dy/2;
    [xd,yd]=meshgrid(xi,yi);
    dxint=5;
    xii=xi(1):dxint:xi(end);
    yii=yi(1):dxint:yi(end);
  [yi1,xi1]=meshgrid(yi,xi);
        [yii,xii]=meshgrid(yii,xii);
        pcaDens60=[];
     pcaDens10=[];
     cn=0;
for objn=1:2
    Ncell=N(objn);
    
    dens=min(10,0.8*Ncell);
   
   
    angx60L23=[];
    angy60L23=[];
    angx60L4=[];
    angy60L4=[];
    angx10L23=[];
    angy10L23=[];
    angx10L4=[];
    angy10L4=[];
    Xaxis=Cs{objn}.Xaxis;
    Yaxis=Cs{objn}.Yaxis;
    Zpk60avg=Cs{objn}.Zpk60avg;
    Zchg60avg=Cs{objn}.Zchg60avg;
    Zden60avg=Cs{objn}.Zden60avg;
    Zlat60avg=Cs{objn}.Zlat60avg;
    Zpk10avg=Cs{objn}.Zpk10avg;
    Zchg10avg=Cs{objn}.Zchg10avg;
    Zden10avg=Cs{objn}.Zden10avg;
    Zlat10avg=Cs{objn}.Zlat10avg;
    ZchgD60avg=Cs{objn}.ZchgD60avg;
    ZdenD60avg=Cs{objn}.ZdenD60avg;
    ZlatD60avg=Cs{objn}.ZlatD60avg;
    Bnd=Cs{objn}.Bnd;
    
    
    
    
     DensMLayer60=cell(Ncell,1);
     ChgMLayer60=cell(Ncell,1);
     PeakMLayer60=cell(Ncell,1);
    LatMlayer60=cell(Ncell,1); 
    DensMLayer10=cell(Ncell,1);
     ChgMLayer10=cell(Ncell,1);
     PeakMLayer10=cell(Ncell,1);
    LatMlayer10=cell(Ncell,1); 
    XaxisLayer=cell(Ncell,1);
       YaxisLayer=cell(Ncell,1);
  k=1;
    for i=1:Ncell

        cn=cn+1;
        indang60=find(Zden60avg{i}>0&Xaxis{i}>=Bnd(i,2)&Xaxis{i}<Bnd(i,3));
        angx60L23=[angx60L23 Xaxis{i}(indang60)];
        angy60L23=[angy60L23 Yaxis{i}(indang60)];
         indang10=find(Zden10avg{i}>0&Xaxis{i}>=Bnd(i,2)&Xaxis{i}<Bnd(i,3));
        angx10L23=[angx10L23 Xaxis{i}(indang10)];
        angy10L23=[angy10L23 Yaxis{i}(indang10)];
        
        indang60=find(Zden60avg{i}>0&Xaxis{i}>=Bnd(i,3)&Xaxis{i}<Bnd(i,4));
        angx60L4=[angx60L4 Xaxis{i}(indang60)];
        angy60L4=[angy60L4 Yaxis{i}(indang60)];
         indang10=find(Zden10avg{i}>0&Xaxis{i}>=Bnd(i,3)&Xaxis{i}<Bnd(i,4));
        angx10L4=[angx10L4 Xaxis{i}(indang10)];
        angy10L4=[angy10L4 Yaxis{i}(indang10)];
        
%         DensV60=[];
%         
%           DensV10=[];
%           
        avg60 = gridavg2(Xaxis{i},Yaxis{i},Zden60avg{i},xtick,ytick,false,Ncell,0);
       
%         for dcircl=0:dx:1500
%          [indx indy]=find(sqrt(xd.^2+yd.^2)<dcircl+dx&sqrt(xd.^2+yd.^2)>=dcircl);
%          vv=[];
%          for iv=1:length(indx)
%          vv=[vv;avg60(indx(iv),indy(iv))];
%          end
%          DensV60=[DensV60;vv];
%     
%         end
        

        avg60=avg60';
         pcaDens60(:,:,cn)=avg60;
        avg2_60=interp2(yi1,xi1,avg60,yii,xii,'linear');
       
%                 currfig= find_figure(sprintf('Singlecelldensity60_%d_%d',objn,k));
%         subplot(3,3,mod(i,9)+1)
%         imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;axis off;
        hold on
%         for ib=1:length(Bnd(1,:))
%             line([-100 100],[Bnd(i,ib) Bnd(i,ib)],'LineWidth',2,'Color',[1 1 1])
%         end
        
        for ib=1:length(Bnd(1,:))-1
            ind{ib}=find(xii(:,1)>=Bnd(i,ib)&xii(:,1)<Bnd(i,ib+1));
            densM60{ib}=avg2_60(ind{ib},:);
            Xmax=max(xii(ind{ib},1));
            Xmin=min(xii(ind{ib},1));
            xcentre=(Xmax+Xmin)/2;
            Xax{ib}=(xii(ind{ib},:)-xcentre)./abs(Bnd(i,ib)-Bnd(i,ib+1)).*mean(abs(Bnd(:,ib)-Bnd(:,ib+1)));
            Yax{ib}=yii(ind{ib},:);
            
        end
        
          %__________________________________________________________________
        
        
        
        avg10 = gridavg2(Xaxis{i},Yaxis{i},Zden10avg{i},xtick,ytick,false,Ncell,0);
        
%           for dcircl=0:dx:1500
%          [indx indy]=find(sqrt(xd.^2+yd.^2)<dcircl+dx&sqrt(xd.^2+yd.^2)>=dcircl);
%           vv=[];
%          for iv=1:length(indx)
%          vv=[vv;avg10(indx(iv),indy(iv))];
%          end
%         
%          DensV10=[DensV10;vv(:)];
%     
%         end
%          pcaDens10=[pcaDens10;DensV10'];

        avg10=avg10';
         pcaDens10(:,:,cn)=avg10;
        avg2_10=interp2(yi1,xi1,avg10,yii,xii,'linear');

        hold on
      
        
        for ib=1:length(Bnd(1,:))-1
            
            densM10{ib}=avg2_10(ind{ib},:);
            
            
        end
        
         C=repmat(avg2_60,[1,1,3]);
      C1=repmat(avg2_60,[1,1,3]);
      C2=repmat(avg2_60,[1,1,3]);
      C2(:,:,1)=avg2_60;
      %[cid1,cid2]=find(avg2_60+avg2_10);
      C2(:,:,2)=avg2_60;
      %C(cid1,cid2,2)=0;
      C2(:,:,3)=0;
      C2(find(C2<0|isnan(C2)|C2>1))=0;
      
      C1(:,:,1)=0;
      %[cid1,cid2]=find(avg2_60+avg2_10);
      C1(:,:,2)=avg2_10;
      %C(cid1,cid2,2)=0;
      C1(:,:,3)=avg2_10;
      C1(find(C1<0|isnan(C1)|C1>1))=0;
      C=(1-C1).*(1-C2);
               
          currfig= find_figure(sprintf('Singlecelldensity_%d_%d',objn,k));
      subplot(3,3,mod(i,9)+1)
        for ib=1:length(Bnd(1,:))
            line([-100 100],[Bnd(i,ib) Bnd(i,ib)],'LineWidth',2,'Color',[0 0 0])
        end
        imagesc(yii(1,:),xii(:,1),C); axis image;axis ij;axis off; 
        
        
        %     colormap(winter)
%         colorbar
        
        
        
        
        
        %     colormap(winter)
%         colorbar
        %------------------------------------------------------------------------
        
        avg60 = gridavg2(Xaxis{i},Yaxis{i},Zpk60avg{i},xtick,ytick,false,Ncell,0);
        
      
        avg60=avg60';
        
        avg2_60=interp2(yi1,xi1,avg60,yii,xii,'linear');
       
%           currfig= find_figure(sprintf('Singlecellpeak60_%d_%d',objn,k));
%         subplot(3,3,mod(i,9)+1)
%         imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij; axis off;
        for ib=1:length(Bnd(1,:))
            line([-100 100],[Bnd(i,ib) Bnd(i,ib)],'LineWidth',2,'Color',[0 0 0])
        end
        %     colormap(winter)
%         colorbar
        
        
        for ib=1:length(Bnd(1,:))-1
            
            PkM60{ib}=avg2_60(ind{ib},:);
            
        end
        
        %------------------------------------------------------------------------
        
        avg10 = gridavg2(Xaxis{i},Yaxis{i},Zpk10avg{i},xtick,ytick,false,Ncell,0);
        
       
        avg10=avg10';
        
        avg2_10=interp2(yi1,xi1,avg10,yii,xii,'linear');
        
%          currfig= find_figure(sprintf('Singlecellpeak10_%d_%d',objn,k));
%         subplot(3,3,mod(i,9)+1)
%         imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;axis off;
       
        %     colormap(winter)
%         colorbar
        
        
        for ib=1:length(Bnd(1,:))-1
            
            PkM10{ib}=avg2_10(ind{ib},:);
            
        end
%           C=repmat(avg2_60,[1,1,3]);
%       C1=repmat(avg2_60,[1,1,3]);
%       C2=repmat(avg2_60,[1,1,3]);
%       C2(:,:,1)=avg2_60;
%       %[cid1,cid2]=find(avg2_60+avg2_10);
%       C2(:,:,2)=avg2_60;
%       %C(cid1,cid2,2)=0;
%       C2(:,:,3)=0;
%       C2(find(C2<0|isnan(C2)|C2>1))=0;
%       
%       C1(:,:,1)=0;
%       %[cid1,cid2]=find(avg2_60+avg2_10);
%       C1(:,:,2)=avg2_10;
%       %C(cid1,cid2,2)=0;
%       C1(:,:,3)=avg2_10;
%       C1(find(C1<0|isnan(C1)|C1>1))=0;
%       C=(1-C1).*(1-C2);
%                
%           currfig= find_figure(sprintf('SinglecellPeak_%d_%d',objn,k));
%       subplot(3,3,mod(i,9)+1)
%  for ib=1:length(Bnd(1,:))
%             line([-100 100],[Bnd(i,ib) Bnd(i,ib)],'LineWidth',2,'Color',[0 0 0])
%         end
%         imagesc(yii(1,:),xii(:,1),C); axis image;axis ij;axis off; 
%         
        %---------------------------------------------------------------
        
        avg60 = gridavg2(Xaxis{i},Yaxis{i},Zchg60avg{i},xtick,ytick,false,Ncell,0);
        

        avg60=avg60';
        
        avg2_60=interp2(yi1,xi1,avg60,yii,xii,'linear');
       
%                currfig= find_figure(sprintf('SinglecellChg60_%d_%d',objn,k));
%         subplot(3,3,mod(i,9)+1) 
%         imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;axis off;
      
        %     colormap(winter)
%         colorbar
        
        
        for ib=1:length(Bnd(1,:))-1
            
            ChgM60{ib}=avg2_60(ind{ib},:);
            
        end
        
        %---------------------------------------------------------------
        
        avg10 = gridavg2(Xaxis{i},Yaxis{i},Zchg10avg{i},xtick,ytick,false,Ncell,0);
        
       
        avg10=avg10';
        
        avg2_10=interp2(yi1,xi1,avg10,yii,xii,'linear');
       
%          currfig= find_figure(sprintf('SinglecellChg10_%d_%d',objn,k));
%         subplot(3,3,mod(i,9)+1)
%         imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;axis off;
%         for ib=1:length(Bnd(1,:))
%             line([-100 100],[Bnd(i,ib) Bnd(i,ib)],'LineWidth',2,'Color',[0 0 0])
%         end
        %     colormap(winter)
%         colorbar
        
        
        for ib=1:length(Bnd(1,:))-1
            
            ChgM10{ib}=avg2_10(ind{ib},:);
            
        end 
%           C=repmat(avg2_60,[1,1,3]);
%       C1=repmat(avg2_60,[1,1,3]);
%       C2=repmat(avg2_60,[1,1,3]);
%       C2(:,:,1)=avg2_60;
%       %[cid1,cid2]=find(avg2_60+avg2_10);
%       C2(:,:,2)=avg2_60;
%       %C(cid1,cid2,2)=0;
%       C2(:,:,3)=0;
%       C2(find(C2<0|isnan(C2)|C2>1))=0;
%       
%       C1(:,:,1)=0;
%       %[cid1,cid2]=find(avg2_60+avg2_10);
%       C1(:,:,2)=avg2_10;
%       %C(cid1,cid2,2)=0;
%       C1(:,:,3)=avg2_10;
%       C1(find(C1<0|isnan(C1)|C1>1))=0;
%       C=(1-C1).*(1-C2);
%                
%           currfig= find_figure(sprintf('SinglecellCharge_%d_%d',objn,k));
%       subplot(3,3,mod(i,9)+1)
%   for ib=1:length(Bnd(1,:))
%             line([-100 100],[Bnd(i,ib) Bnd(i,ib)],'LineWidth',2,'Color',[0 0 0])
%         end
%         imagesc(yii(1,:),xii(:,1),C); axis image;axis ij;axis off; 
        
        
        
        %---------------------------------------------------------------
        avg60 = gridavg2(Xaxis{i},Yaxis{i},Zlat60avg{i},xtick,ytick,false,Ncell,0);
        

        avg60=avg60';
        
        avg2_60=interp2(yi1,xi1,avg60,yii,xii,'linear');
       
%                 currfig= find_figure(sprintf('SinglecellLatency60_%d_%d',objn,k));
%        subplot(3,3,mod(i,9)+1)
%         imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;axis off;
%         for ib=1:length(Bnd(1,:))
%             line([-100 100],[Bnd(i,ib) Bnd(i,ib)],'LineWidth',2,'Color',[0 0 0])
%         end
        %     colormap(winter)
%         colorbar
        
        
        for ib=1:length(Bnd(1,:))-1
            
            LatM60{ib}=avg2_60(ind{ib},:);
            
        end
       
       
       
        %---------------------------------------------------------------
        avg10 = gridavg2(Xaxis{i},Yaxis{i},Zlat10avg{i},xtick,ytick,false,Ncell,0);
        
        
        avg10=avg10';
        
        avg2_10=interp2(yi1,xi1,avg10,yii,xii,'linear');
     
%         currfig= find_figure(sprintf('SinglecellLatency10_%d_%d',objn,k));
%         subplot(3,3,mod(i,9)+1)
%         imagesc(yii(1,:),xii(:,1),avg2); axis image;axis ij;axis off;
        %     colormap(winter)
%         colorbar
        
        for ib=1:length(Bnd(1,:))-1
            
            LatM10{ib}=avg2_10(ind{ib},:);
            
        end
        
%          C=repmat(avg2_60,[1,1,3]);
%       C1=repmat(avg2_60,[1,1,3]);
%       C2=repmat(avg2_60,[1,1,3]);
%       C2(:,:,1)=avg2_60;
%       %[cid1,cid2]=find(avg2_60+avg2_10);
%       C2(:,:,2)=avg2_60;
%       %C(cid1,cid2,2)=0;
%       C2(:,:,3)=0;
%       C2(find(C2<0|isnan(C2)|C2>1))=0;
%       
%       C1(:,:,1)=0;
%       %[cid1,cid2]=find(avg2_60+avg2_10);
%       C1(:,:,2)=avg2_10;
%       %C(cid1,cid2,2)=0;
%       C1(:,:,3)=avg2_10;
%       C1(find(C1<0|isnan(C1)|C1>1))=0;
%       C=(1-C1).*(1-C2);
%                
%           currfig= find_figure(sprintf('SinglecellLatency_%d_%d',objn,k));
%       subplot(3,3,mod(i,9)+1)
%         for ib=1:length(Bnd(1,:))
%             line([-100 100],[Bnd(i,ib) Bnd(i,ib)],'LineWidth',2,'Color',[0 0 0])
%         end

%         imagesc(yii(1,:),xii(:,1),C); axis image;axis ij;axis off; 
        %__________________________________________________________________
        
        
        DensMLayer60{i}=densM60;
        PeakMLayer60{i}=PkM60;
        ChgMLayer60{i}=ChgM60;
        LatMlayer60{i}=LatM60;
         DensMLayer10{i}=densM10;
        PeakMLayer10{i}=PkM10;
        ChgMLayer10{i}=ChgM10;
        LatMlayer10{i}=LatM10;
        XaxisLayer{i}=Xax;
        YaxisLayer{i}=Yax;
        if mod(i,9)==0
            k=k+1;
        end

    end
    LayerAvg{objn}.DensMLayer60=DensMLayer60;
    LayerAvg{objn}.PeakMLayer60=PeakMLayer60; 
     LayerAvg{objn}.ChgMLayer60=ChgMLayer60;
    LayerAvg{objn}.LatMlayer60=LatMlayer60; 
     LayerAvg{objn}.DensMLayer10=DensMLayer10;
    LayerAvg{objn}.PeakMLayer10= PeakMLayer10; 
     LayerAvg{objn}.ChgMLayer10=ChgMLayer10;
    LayerAvg{objn}.LatMlayer10=LatMlayer10;
      LayerAvg{objn}.XaxisLayer=XaxisLayer;
    LayerAvg{objn}.YaxisLayer=YaxisLayer;
    LayerAvg{objn}.angx60L23=angx60L23;
    LayerAvg{objn}.angx10L23=angx10L23;
    LayerAvg{objn}.angy60L23=angy60L23;
    LayerAvg{objn}.angy10L23=angy10L23;
     LayerAvg{objn}.angx60L4=angx60L4;
    LayerAvg{objn}.angx10L4=angx10L4;
    LayerAvg{objn}.angy60L4=angy60L4;
    LayerAvg{objn}.angy10L4=angy10L4;
    
end

for objn=1:2;
Bnd=Cs{objn}.Bnd;
Ncell=N(objn);
density=0.1; 
angx60L23=LayerAvg{objn}.angx60L23;
angy60L23=LayerAvg{objn}.angy60L23;
angx10L23=LayerAvg{objn}.angx10L23;
angy10L23=LayerAvg{objn}.angy10L23;
ang60L23=atan(angy60L23./abs(angx60L23))+pi./2;
ang10L23=atan(angy10L23./abs(angx10L23))+pi./2;

angx60L4=LayerAvg{objn}.angx60L4;
angy60L4=LayerAvg{objn}.angy60L4;
angx10L4=LayerAvg{objn}.angx10L4;
angy10L4=LayerAvg{objn}.angy10L4;
ang60L4=atan(angy60L4./abs(angx60L4))+pi./2;
ang10L4=atan(angy10L4./abs(angx10L4))+pi./2;


theta1=-0.5:0.1:0.5;
ahist60L23=hist(angy60L23./1000,theta1)./length(angy60L23);
ahist60L23_1=zeros(size(ahist60L23));
ahist60L23_1(find(ahist60L23>0.15))=1;
ahist60L4L=hist(angy60L4./1000,theta1)./length(angy60L4);
ahist60L4L_1=zeros(size(ahist60L4L));
ahist60L4L_1(find(ahist60L4L>0.15))=1;

r1=(0:1:1)';
[X1,Y1]=meshgrid(theta1,r1)
C = meshgrid(ahist60L23,r1);
find_figure(sprintf('AngledensityDistribution60%i',objn))
subplot(2,1,1)
pcolor(X1,Y1,C)
subplot(2,1,2)



theta = pi*(0:2:18)/18;
ahist60L4=hist(ang60L4,theta)./length(ang60L4);
ahist60L4_1=zeros(size(ahist60L4));
ahist60L4_1(find(ahist60L4>0.15))=1;
n = 6;
r = (3:n)'/n;

X = r*cos(theta);
Y = r*sin(theta);
C = meshgrid(ahist60L4,r);
pcolor(X,Y,C)
axis equal tight 

C = meshgrid(ahist60L23,r1);
find_figure(sprintf('LinedensityDistribution60%i',objn))
subplot(2,1,1)
pcolor(X1,Y1,C)
subplot(2,1,2)
C = meshgrid(ahist60L4L,r1);
pcolor(X1,Y1,C)

find_figure(sprintf('AngledensityDistribution60%i_cutoff',objn))
subplot(2,1,1)
C = meshgrid(ahist60L23_1,r1);
pcolor(X1,Y1,C)
subplot(2,1,2)


C = meshgrid(ahist60L4_1,r);
pcolor(X,Y,C)
axis equal tight 

% theta=-0.5:0.1:0.5;
ahist10L23=hist(angy10L23./1000,theta1)./length(angy10L23);
ahist10L23_1=zeros(size(ahist10L23));
ahist10L23_1(find(ahist10L23>0.15))=1;


ahist10L4L=hist(angy10L4./1000,theta1)./length(angy10L4);
ahist10L4L_1=zeros(size(ahist10L4L));
ahist10L4L_1(find(ahist10L4L>0.15))=1;
% r1=(0:1:1)';
% [X,Y]=meshgrid(theta,r)
C = meshgrid(ahist10L23,r1);
find_figure(sprintf('AngledensityDistribution10%i',objn))
subplot(2,1,1)
pcolor(X1,Y1,C)
subplot(2,1,2)

% theta = pi*(0:2:18)/18;
ahist10L4=hist(ang10L4,theta)./length(ang10L4);
ahist10L4_1=zeros(size(ahist10L4));
ahist10L4_1(find(ahist10L4>0.15))=1;
% n = 6;
% r = (3:n)'/n;
% 
% X = r*cos(theta);
% Y = r*sin(theta);
C = meshgrid(ahist10L4,r);
pcolor(X,Y,C)
axis equal tight 


C = meshgrid(ahist10L23,r1);
find_figure(sprintf('LinedensityDistribution10%i',objn))
subplot(2,1,1)
pcolor(X1,Y1,C)
subplot(2,1,2)
C = meshgrid(ahist10L4L,r1);
pcolor(X1,Y1,C)

C = meshgrid(ahist10L23_1,r1);
find_figure(sprintf('AngledensityDistribution10%icutoff',objn))
subplot(2,1,1)
pcolor(X1,Y1,C)
subplot(2,1,2)

% theta = pi*(0:2:18)/18;
% ahist10L4=hist(ang10L4,theta)./length(ang10L4);
% ahist10L4_1=ahist10L4>0.15;
% n = 6;
% r = (3:n)'/n;
% 
% X = r*cos(theta);
% Y = r*sin(theta);
C = meshgrid(ahist10L4_1,r);
pcolor(X,Y,C)
axis equal tight 
end
%avg Layer1
% for j=1:length(Bnd(1,:))-1
%    Zdens60=[];
% Zpk60=[];
% Zchg60=[];
% Zdens10=[];
% Zpk10=[];
% Zchg10=[];
% X=[];
% Y=[]; 
% xtick=LayerAvg{objn}.YaxisLayer{1}{1}(1,:)-dxint/2;
% ytick=-1*mean(abs(Bnd(:,j)-Bnd(:,j+1)))/2-dxint/2:dxint:mean(abs(Bnd(:,j)-Bnd(:,j+1)))/2-dxint/2;
% xi=xtick+0.5;
% yi=ytick+0.5;
% for i=1:Ncell
%   X=[X;LayerAvg{objn}.XaxisLayer{i}{j}(:)];
%   Y=[Y;LayerAvg{objn}.YaxisLayer{i}{j}(:)];
%   Zdens60=[Zdens60;LayerAvg{objn}.DensMLayer60{i}{j}(:)];
%   Zchg60=[Zchg60;LayerAvg{objn}.ChgMLayer60{i}{j}(:)];
%   Zpk60=[Zpk60;LayerAvg{objn}.PeakMLayer60{i}{j}(:)];
%   Zdens10=[Zdens10;LayerAvg{objn}.DensMLayer10{i}{j}(:)];
%   Zchg10=[Zchg10;LayerAvg{objn}.ChgMLayer10{i}{j}(:)];
%   Zpk10=[Zpk10;LayerAvg{objn}.PeakMLayer10{i}{j}(:)];
% end
% 
% avg = gridavg2(Y,X,Zdens60,xtick,ytick,false,Ncell,0);
% Ind60=find(avg<density);
% 
% currfig= find_figure(sprintf('AvgDensity60_object%b',objn));
% subplot(length(Bnd(1,:))-1,1,j)
% imagesc(xi,yi,avg); axis image;axis ij;
% currfig= find_figure(sprintf('Avgchg60_object%b',objn));
% avg = gridavg2(Y,X,Zchg60,xtick,ytick,false,Ncell,0);
% avg(Ind60)=0;
% subplot(length(Bnd(1,:))-1,1,j)
% imagesc(xi,yi,avg); axis image;axis ij;
% 
% currfig= find_figure(sprintf('Avgpk60_object%b',objn));
% avg = gridavg2(Y,X,Zpk60,xtick,ytick,false,Ncell,0);
% avg(Ind60)=0;
% subplot(length(Bnd(1,:))-1,1,j)
% imagesc(xi,yi,avg); axis image;axis ij;
% 
% 
% avg = gridavg2(Y,X,Zdens10,xtick,ytick,false,Ncell,0);
% Ind10=find(avg<density);
% currfig= find_figure(sprintf('AvgDensity10_object%b',objn));
% subplot(length(Bnd(1,:))-1,1,j)
% imagesc(xi,yi,avg); axis image;axis ij;
% currfig= find_figure(sprintf('Avgchg10_object%b',objn));
% avg = gridavg2(Y,X,Zchg10,xtick,ytick,false,Ncell,0);
% avg(Ind10)=0;
% subplot(length(Bnd(1,:))-1,1,j)
% imagesc(xi,yi,avg); axis image;axis ij;
% 
% currfig= find_figure(sprintf('Avgpk10_object%b',objn));
% avg = gridavg2(Y,X,Zpk10,xtick,ytick,false,Ncell,0);
% avg(Ind10)=0;
% subplot(length(Bnd(1,:))-1,1,j)
% imagesc(xi,yi,avg); axis image;axis ij;
% 
% end



    
    
    
