folderMeanpath='/Users/lindameng/Documents/study/kawens/Norm/cellattach/eventwindow50';
a=dir(fullfile(folderMeanpath,'*.mat'));
Cellfile=size(a);
eventwindow=50;
Dist80=[];
Dist80Event=[];
for j= 1:1:Cellfile
    Cg= load(fullfile(folderMeanpath,a(j).name));
    if isfield(Cg,'cells')
    cells=Cg.cells;
    header=cells.header;
     dx = cells.header.PatternSpacing(1);
    xdist=0:dx:1000;
    sel=cells.Latency>0;
    ind=find(sel);
        xaxis=cells.header.StimCoordinates(1,ind);
    yaxis=cells.header.StimCoordinates(2,ind);
     Totalevent=sum(sel);
    Totalevent60=sum(cells.Latency>0&cells.Latency<=eventwindow);
    ind60=find(cells.Latency>0&cells.Latency<=eventwindow);
    xaxis60=cells.header.StimCoordinates(1,ind60);
    yaxis60=cells.header.StimCoordinates(2,ind60);
%         xsoma=cells.header.Soma1Coordinates(1,1);
%     ysoma=cells.header.Soma1Coordinates(1,2);
    if length(header.Soma1Coordinates)
    xsoma=header.Soma1Coordinates(1,1);
    ysoma=header.Soma1Coordinates(1,2);


else
   
   h=figure(300)
    [imd imdIdx] = sort(header.ImageData(:));
    imd(imdIdx) = floor(linspace(0,256-eps,numel(header.ImageData)));
    imd = uint8(reshape(imd, size(header.ImageData)));
    
    
    ImageScale=header.ImageScale;
    
    
    image([-1 1].*ImageScale(1)/2, [-1 1].*ImageScale(2)/2, repmat(imd,[1 1 3]));
    axis xy;
    set(gca,'XTick',[-1 -.5 0 .5 1].*ImageScale(1)/2, 'YTick',[-1 -.5 0 .5 1].*ImageScale(2)/2, ...
        'DataAspectRatio',[ 1 1 1]);
    [x,y]=ginput(1);
    header.Soma1Coordinates=zeros(1,2);
    header.Soma1Coordinates(1,1)=x;
    header.Soma1Coordinates(1,2)=y;
    xsoma=x;
    ysoma=y;
  
    close(h)
end
    
    dist=sqrt((xaxis-xsoma).^2+(yaxis-ysoma).^2);
    
     dist60=sqrt((xaxis60-xsoma).^2+(yaxis60-ysoma).^2);
       perc=zeros(length(xdist),1);
    perc60=zeros(length(xdist),1);
    for i=1:1:length(xdist)
        perc(i)=sum(dist<=xdist(i))/Totalevent*100;
        perc60(i)=sum(dist60<=xdist(i))/Totalevent60*100;
        
    end
    cells.perc=perc;
    cells.perc60=perc60;
    cells.xdist=xdist;
    ii=2:1:length(xdist);
    ind80=find(perc(ii-1)<80&perc(ii)>=80);
    ind80E=find(perc60(ii-1)<80&perc60(ii)>=80);
    if ~isempty(ind80)
    Dist80=[Dist80;xdist(ind80(1))];
    
    else
        Dist80=[Dist80;0];
       
    end
        if ~isempty(ind80E)
    
    Dist80Event=[Dist80Event;xdist(ind80E(1))];
    else
     
        Dist80Event=[Dist80Event;0];
    end
    datafolder=fullfile(folderMeanpath,a(j).name);
    save(datafolder,'cells');
    end
end
fileDist80=fullfile(folderMeanpath,'Dist80Perc');
save(fileDist80,'Dist80','Dist80Event');
hist(Dist80)
