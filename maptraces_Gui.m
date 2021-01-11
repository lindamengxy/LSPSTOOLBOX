function handles= maptraces_Gui(data,header,mapname,foname,Tp,winwidth,axestouse,events)
%filefolder is the folder name where your data is. exp: ''/Users/lindameng/Documents/study/meng/012213/bs0001''
%  mapname is the data's name without .dat or .mat. exp: 'Map_130122_1437',

% nPre = 4;
% dN = 50;
% offset = 1001-dN*nPre;
% N = min((size(data,1)-offset+1)/dN, 100+nPre);
% N1=floor(winwidth*10/dN);
%N = min((size(data_,1)-offset+1)/dN, 100+nPre);


%meds = median(data_(1:800, :));
%data_ = squeeze(sum(reshape(data_(offset:offset+N*dN-1,:)-repmat(meds,[N*dN 1]),[dN N size(data,2)])/1e4));
% meds = median(data(1:800, :));
% data_ = squeeze(sum(reshape(data(offset:offset+N*dN-1,:)-repmat(meds,[N*dN 1]),[dN N size(data,2)])/1e4));
% data_=data_(1:N1+nPre,:);
% if header.holdingPotential<0
%     mData = min(data_);
% else
%     mData = max(data_);
% end
% 
% 
% [pixVals pixValIndices] = sort(data_(:));
% N=size(data_,1);
% if header.holdingPotential<0
%     cLims = pixVals( [10*N round(.3*length(pixVals))]);
% else
%     cLims = pixVals( [round(.7*length(pixVals)) length(pixVals)-10*N+1 ]);
% end
% 
% mpv = mean(pixVals(round(.2*length(pixVals)):round(.8*length(pixVals))));
% spv = std(pixVals(round(.2*length(pixVals)):round(.8*length(pixVals))));
% z = (pixVals - mpv)./spv;
% 
% excBounds = [1 max(find(z<-5))];
% inhBounds = [min(find(z>5)) length(z)];
% 
% %define direct by amplituder criterion
% directmultiplier=6;%4   %used 6 for most data
% [directVal directBound] = min(abs(pixVals-directmultiplier*pixVals(round(mean(excBounds)))));
% [directinhVal directinhBound] = min(abs(pixVals-directmultiplier*pixVals(round(mean(inhBounds)))));
% 
% %cLims1 = pixVals([min(find(z>-100)) max(find(z<-3))]);
%cLims2 = pixVals([min(find(z>3)) max(find(z))]);

% if header.holdingPotential<0
%     c = pixVals([directBound excBounds(2)]);
% else
%     cLims = pixVals([inhBounds(1) directinhBound]);
% end
% directBound = pixVals(directBound);
% directBound=-1000;   %dont worry aboiut amplitude! just basic time! 1/21/11 PK
%
% directIndices = (min(data_) < directBound);
% mDataVals = (mData-cLims(1))./(cLims(2)-cLims(1));
% mDataVals(mDataVals > 1) = 1;
% mDataVals(mDataVals < 0) = 0; %2
% mDataVals = floor(mDataVals*256)+1;
% mDataVals(mDataVals==257) = 256;
% mDataVals(directIndices) = 257; % > 257
[imd imdIdx] = sort(header.ImageData(:));
imd(imdIdx) = floor(linspace(0,256-eps,numel(header.ImageData)));
imd = uint8(reshape(imd, size(header.ImageData)));
stim_start=find(header.stimulatorOutput(4,:));


rsrate=1000; %resampling  sample rate is 10000  => 5x undersampled
resampdata = resample(data,rsrate, size(data,1));
resampdata1 = resample(data,rsrate, size(data,1));
resampdata = resampdata(10:end,:);
if header.holdingPotential<0
resampdata = -(resampdata - repmat(max(resampdata),[rsrate-9 1]))./(max(max(resampdata))-min(min(resampdata)));
else
 resampdata = -(resampdata - repmat(min(resampdata),[rsrate-9 1]))./(max(max(resampdata))-min(min(resampdata)));
end
[pts,stims]=size(data);

 set(gcf,'CurrentAxes',axestouse);
    cla
    ImageScale=header.ImageScale;
    image([-1 1].*ImageScale(2)/2, [-1 1].*ImageScale(1)/2, repmat(imd',[1 1 3]));
    % axis image;
    set(gca,'XTick',[-1 -.5 0 .5 1].*ImageScale(2)/2, 'YTick',[-1 -.5 0 .5 1].*ImageScale(1)/2, ...
        'DataAspectRatio',[ 1 1 1]);
    hold on;
    %     cmap = jet(256);
    %     cmap(257,:) = [0 0 0];
    xRatio = 1;
    yRatio = 1;
    %dX = header.PatternSpacing(1);
    %dY = header.PatternSpacing(2);
    %flipped x and y for paper figure 1/2011 PK
    dX = header.PatternSpacing(2);
    dY = header.PatternSpacing(1);
    t_ = linspace(0, 1, size(resampdata,1));
    
    for k = 1:size(resampdata,2)
        %rects(k) = rectangle('Position', [header.StimCoordinates(:,k)'-[xRatio.*dX yRatio.*dY]/2 [xRatio.*dX yRatio.*dY]],'EdgeColor',[0 0 0]);
        %plots(k) = plot(t_ .* xRatio .* dX + header.StimCoordinates(1,k), (resampdata(:,k)-.5) .* yRatio * dY + header.StimCoordinates(2,k), '-', ...
        % 'Color', cmap(mDataVals(k),:));
        %plots(k) = plot(t_ .* xRatio .* dX + header.StimCoordinates(1,k), smooth((resampdata(:,k)-.5),11) .* yRatio * dY + header.StimCoordinates(2,k), '-', ...
        %   'Color', cmap(mDataVals(k),:));
        
        fevent=find(events.flag{k}(:)>0);
        if length(fevent)>0
            if events.flag{k}(fevent(1))==1
                plots(k) = plot(t_(1:500) .* xRatio .* dY + header.StimCoordinates(2,k),smooth((resampdata(1:500,k)),11) .* yRatio * dX + header.StimCoordinates(1,k), 'k-');
            elseif events.flag{k}(fevent(1))==2
                plots(k) = plot(t_(1:500).* xRatio .* dY + header.StimCoordinates(2,k),smooth((resampdata(1:500,k)),11) .* yRatio * dX + header.StimCoordinates(1,k), 'r-');
            end
        else         plots(k) = plot(t_(1:500) .* xRatio .* dY + header.StimCoordinates(2,k),smooth((resampdata(1:500,k)),11) .* yRatio * dX + header.StimCoordinates(1,k), 'b-');
        end
        ImageScale=header.ImageScale;
        
        set(gca,'XTick',[-1 -.5 0 .5 1].*ImageScale(2)/2, 'YTick',[-1 -.5 0 .5 1].*ImageScale(1)/2, ...
            'DataAspectRatio',[ 1 1 1]);
        
    end
    fname = regexp(header.dataFilename,'\\Map_(?<fname>.+)\.dat','tokens');
    fnamestr=char(fname{1});
    title(sprintf('%s\n%s %imV ','trace',fnamestr,header.holdingPotential));
    hold off;

        foldername=fullfile(foname);
        if (exist(foldername) == 0)
            mkdir (foldername);
        end
%         f1=fullfile(foldername,experiment);
%         if (exist(f1) == 0)
%             mkdir (f1);
%         end
%         f2=fullfile(f1,cell);
%         if (exist(f2) == 0)
%             mkdir (f2);
%         end
        if Tp==2
            typefolder='PTX';
        elseif Tp==3
            typefolder='APV_PTX';
        elseif Tp==1
            typefolder='highMg';
        elseif Tp==4
            typefolder='TTX';
        else
            typefolder=handles.additiondrugs;
        end
        f4=fullfile(foldername,typefolder);
        if (exist(f4) == 0)
            mkdir (f4);
        end

filesname1=fullfile(f4,sprintf('Maptrace_%s_eventwindow%i.eps',mapname,winwidth));
filesname2=fullfile(f4,sprintf('Maptrace_%s_eventwindow%i.fig',mapname,winwidth));
filesname3=fullfile(f4,sprintf('Maptrace_%s_eventwindow%i.png',mapname,winwidth));


%guidata(hb, handles);

saveImages=1;
if saveImages == 1
    
    
  currfig=find_figure('Trace_DIC');
    clf;
    ImageScale=header.ImageScale;
    image([-1 1].*ImageScale(2)/2, [-1 1].*ImageScale(1)/2, repmat(imd',[1 1 3]));
     axis image;
    set(gca,'XTick',[-1 -.5 0 .5 1].*ImageScale(2)/2, 'YTick',[-1 -.5 0 .5 1].*ImageScale(1)/2, ...
        'DataAspectRatio',[ 1 1 1]);
    hold on;
    %     cmap = jet(256);
    %     cmap(257,:) = [0 0 0];
    xRatio = 1;
    yRatio = 1;
    %dX = header.PatternSpacing(1);
    %dY = header.PatternSpacing(2);
    %flipped x and y for paper figure 1/2011 PK
    dX = header.PatternSpacing(2);
    dY = header.PatternSpacing(1);
    t_ = linspace(0, 1, size(resampdata,1));
    
    for k = 1:size(resampdata,2)
        %rects(k) = rectangle('Position', [header.StimCoordinates(:,k)'-[xRatio.*dX yRatio.*dY]/2 [xRatio.*dX yRatio.*dY]],'EdgeColor',[0 0 0]);
        %plots(k) = plot(t_ .* xRatio .* dX + header.StimCoordinates(1,k), (resampdata(:,k)-.5) .* yRatio * dY + header.StimCoordinates(2,k), '-', ...
        % 'Color', cmap(mDataVals(k),:));
        %plots(k) = plot(t_ .* xRatio .* dX + header.StimCoordinates(1,k), smooth((resampdata(:,k)-.5),11) .* yRatio * dY + header.StimCoordinates(2,k), '-', ...
        %   'Color', cmap(mDataVals(k),:));
        
        fevent=find(events.flag{k}(:)>0);
        if length(fevent)>0
            if events.flag{k}(fevent(1))==1
                plots(k) = plot(t_(1:500) .* xRatio .* dY + header.StimCoordinates(2,k),smooth((resampdata(1:500,k)),11) .* yRatio * dX + header.StimCoordinates(1,k), 'k-');
            elseif events.flag{k}(fevent(1))==2
                plots(k) = plot(t_(1:500).* xRatio .* dY + header.StimCoordinates(2,k),smooth((resampdata(1:500,k)),11) .* yRatio * dX + header.StimCoordinates(1,k), 'r-');
            end
        else         plots(k) = plot(t_(1:500) .* xRatio .* dY + header.StimCoordinates(2,k),smooth((resampdata(1:500,k)),11) .* yRatio * dX + header.StimCoordinates(1,k), 'b-');
        end
        ImageScale=header.ImageScale;
        
        set(gca,'XTick',[-1 -.5 0 .5 1].*ImageScale(1)/2, 'YTick',[-1 -.5 0 .5 1].*ImageScale(2)/2, ...
            'DataAspectRatio',[ 1 1 1]);
        
    end
    fname = regexp(header.dataFilename,'\\Map_(?<fname>.+)\.dat','tokens');
    fnamestr=char(fname{1});
    title(sprintf('%s\n%s %imV ','trace',fnamestr,header.holdingPotential));
    hold off;
    
    print(currfig,'-depsc2',filesname1);
    saveas(currfig,filesname2,'fig');
    print(currfig,'-dpng',filesname3);
    close(currfig);
end



end

