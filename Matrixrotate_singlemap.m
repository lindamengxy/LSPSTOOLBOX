function [stimcoordinates header]=Matrixrotate_singlemap(header,flipimg,flipimg2)
%shift and move the matrix to make sure the cell that patched is in (0,0)
rotate_angle=-header.SpatialRotation/180*3.14;

if length(header.Soma1Coordinates)
    somax=header.Soma1Coordinates(1,1);
    somay=header.Soma1Coordinates(1,2);


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
    title('click the patched cell')
    [x,y]=ginput(1);
    header.Soma1Coordinates=zeros(1,2);
    header.Soma1Coordinates(1,1)=x;
    header.Soma1Coordinates(1,2)=y;
    somax=x;
    somay=y;
    close(h)
end
% somax=header.Soma1Coordinates(1,1);
% somay=header.Soma1Coordinates(1,2);

stimcoordinates=header.StimCoordinates;
%shift cell to (0,0)
soma=[somax;somay];
stimcoordinates=stimcoordinates-repmat(soma,1,size(stimcoordinates,2));
%rotate
rotatmatrix=[cos(rotate_angle),-sin(rotate_angle);sin(rotate_angle), cos(rotate_angle)];


stimcoordinates=rotatmatrix*stimcoordinates;

if flipimg
   stimcoordinates(2,:)=-stimcoordinates(2,:);
end
if flipimg2
   stimcoordinates(1,:)=-stimcoordinates(1,:);
end

end

