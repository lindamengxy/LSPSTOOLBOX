folderMeanpath='/Users/lindameng/Documents/study/meng_code2/silent/meanValueeachcell_eventwindow50_direct8';
a=dir(fullfile(folderMeanpath,'*.mat'));
Cellfile=size(a);
Number_cell=0;
ChgBeforeAPV60=[];
ChgAfterAPV60=[];
ChgBeforeAPV50=[];
ChgAfterAPV50=[];
DistBeforeAPV60=[];
DistAfterAPV60=[];
DistBeforeAPV50=[];
DistAfterAPV50=[];
for i= 1:1:Cellfile
    Cg= load(fullfile(folderMeanpath,a(i).name));
    if isfield(Cg,'Cellrecord')
        Cellrecord=Cg.Cellrecord;
        for fn=1:Cg.Num_cell
             if isfield(Cellrecord{fn},'meanareaapvInd50')||isfield(Cellrecord{fn},'meanareaapvInd60')
                 Number_cell=Number_cell+1
                 if isfield(Cellrecord{fn},'meanareaapvInd60')
                     x0=Cellrecord{fn}.SomaCoordinates(1,1);
                     y0=Cellrecord{fn}.SomaCoordinates(1,2);
                   Ind60=find(Cellrecord{fn}.meanflagapvInd60>0);
                   ChgAfterAPV60=[ChgAfterAPV60,Cellrecord{fn}.meanareaapvInd60(Ind60)'];
                   x=Cellrecord{fn}.StimCoordinates(1,Ind60);
                   y=Cellrecord{fn}.StimCoordinates(2,Ind60);
                   DistAfterAPV60=[DistAfterAPV60,sqrt((x-x0).^2+(y-y0).^2)];
                 end
                     if isfield(Cellrecord{fn},'meanareaapvInd50')
                          x0=Cellrecord{fn}.SomaCoordinates(1,1);
                     y0=Cellrecord{fn}.SomaCoordinates(1,2);
                   Ind50=find(Cellrecord{fn}.meanflagapvInd50>0);
                   ChgAfterAPV50=[ChgAfterAPV50,Cellrecord{fn}.meanareaapvInd50(Ind50)'];
                   x=Cellrecord{fn}.StimCoordinates(1,Ind50);
                   y=Cellrecord{fn}.StimCoordinates(2,Ind50);
                   DistAfterAPV50=[DistAfterAPV50,sqrt((x-x0).^2+(y-y0).^2)];
                     end
                 if isfield(Cellrecord{fn},'meanareaPtxInd60')
                      x0=Cellrecord{fn}.SomaCoordinates(1,1);
                     y0=Cellrecord{fn}.SomaCoordinates(1,2);
                   Ind60=find(Cellrecord{fn}.meanflagPtxInd60>0);
                   ChgBeforeAPV60=[ChgBeforeAPV60,Cellrecord{fn}.meanareaPtxInd60(Ind60)'];
                   x=Cellrecord{fn}.StimCoordinates(1,Ind60);
                   y=Cellrecord{fn}.StimCoordinates(2,Ind60);
                   DistBeforeAPV60=[DistBeforeAPV60,sqrt((x-x0).^2+(y-y0).^2)];
                 end
                     if isfield(Cellrecord{fn},'meanareaPtxInd50')
                          x0=Cellrecord{fn}.SomaCoordinates(1,1);
                     y0=Cellrecord{fn}.SomaCoordinates(1,2);
                   Ind50=find(Cellrecord{fn}.meanflagPtxInd50>0);
                   ChgBeforeAPV50=[ChgBeforeAPV50,Cellrecord{fn}.meanareaPtxInd50(Ind50)'];
                   x=Cellrecord{fn}.StimCoordinates(1,Ind50);
                   y=Cellrecord{fn}.StimCoordinates(2,Ind50);
                   DistBeforeAPV50=[DistBeforeAPV50,sqrt((x-x0).^2+(y-y0).^2)];
                     end
             end
        end
    end
end
    figure(1)
    scatter(log(DistAfterAPV60),log(ChgAfterAPV60))
    figure(2)
      scatter(log(DistAfterAPV60),log(ChgAfterAPV60))
      figure(3)
    scatter(log(DistBeforeAPV60),log(ChgBeforeAPV60))
    figure(4)
      scatter(log(DistBeforeAPV50),log(ChgBeforeAPV50))
      
      

                     