function [filters openFilter] = generateFilters(segPara,doDraw)

Nsize  = 5;
Nangle = 3;

Nf = Nsize + Nangle*Nsize;

filters = zeros(50,50,Nf);
[x y] = meshgrid(1:50,1:50);
d=sqrt( (x-25).^2 + (y-25).^2 );

r = linspace(segPara.minSize,segPara.maxSize,Nsize);

%circular filters
for i=1:Nsize
       
   circle=zeros(50,50);
   
   circle(d<r(i)+5)= -1/ ( (r(i)+5)^2 - (r(i)+3)^2 ); 
   circle(d<r(i)+3)= -1/ ( (r(i)+3)^2 - (r(i)+1)^2 );  
   circle(d<r(i)+1)= -2/ ( (r(i)+3)^2 - (r(i)+1)^2 ); 
   circle(d<r(i))= 3/r(i)^2; 
        
   filters(:,:,i) = circle; 

   
   if doDraw
       clf
       imagesc(circle)
       drawnow
   end
   
end

%% elliptic filters

th = 360/Nangle:360/Nangle:360;

for k=1:Nangle
    for i=1:Nsize

       circle=zeros(50,50);

       circle(d<r(i)+5)= -1/ ( (r(i)+5)^2 - (r(i)+3)^2 ); 
       circle(d<r(i)+3)= -1/ ( (r(i)+3)^2 - (r(i)+1)^2 );  
       circle(d<r(i)+1)= -2/ ( (r(i)+3)^2 - (r(i)+1)^2 ); 
       circle(d<r(i))= 3/r(i)^2; 
       
       circle = circle(11:39,:);


       tform = maketform('affine',[1 0 0; 0 1 0; 0 0 1]);
       ell = imtransform(circle,tform,'Size',[50 50]);
       ell = imrotate(ell,th(k),'crop');

       filters(:,:,i+Nsize*k) = ell;

       
       
       if doDraw
       imagesc(ell)
       drawnow
       end

    end
end


%%


%for open
circle=zeros(50,50);
circle(d<segPara.openSize)= 1; 

openFilter = circle;
