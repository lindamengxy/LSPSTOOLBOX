function Boundry=boundary(header,flipimg,BD)
%shift and move the matrix to make sure the cell that patched is in (0,0)
rotate_angle=-header.SpatialRotation/180*3.14;
currfig=find_figure('Boundary')
clf;
    [imd imdIdx] = sort(header.ImageData(:));
    imd(imdIdx) = floor(linspace(0,256-eps,numel(header.ImageData)));
    imd = uint8(reshape(imd, size(header.ImageData)));
    
    
    ImageScale=header.ImageScale;
    
    
    image([-1 1].*ImageScale(1)/2, [-1 1].*ImageScale(2)/2, repmat(imd,[1 1 3]));
    hold on
    axis xy;
    set(gca,'XTick',[-1 -.5 0 .5 1].*ImageScale(1)/2, 'YTick',[-1 -.5 0 .5 1].*ImageScale(2)/2, ...
        'DataAspectRatio',[ 1 1 1]);
    title('click five points for the boudary [pia L2/3->L4 L4->L5/6 L5/6->subplate whitematter]') 
    [bdx bdy]=ginput(BD);
    plot(bdx,bdy,'r*')
   
if length(header.Soma1Coordinates)
    somax=header.Soma1Coordinates(1,1);
    somay=header.Soma1Coordinates(1,2);


else
    title('click the cell you patch')
    hold off
    [x,y]=ginput(1);
    header.Soma1Coordinates=zeros(1,2);
    header.Soma1Coordinates(1,1)=x;
    header.Soma1Coordinates(1,2)=y;
    somax=x;
    somay=y;
    
end
close(currfig)
%shift cell to (0,0)
soma=[somax;somay];
Boundry=[bdx bdy]';

Boundry=Boundry-repmat(soma,1,size(Boundry,2));
%rotate
rotatmatrix=[cos(rotate_angle),-sin(rotate_angle);sin(rotate_angle), cos(rotate_angle)];

Boundry=rotatmatrix*Boundry;

if flipimg
  Boundry(2,:)=-Boundry(2,:);
end

end

