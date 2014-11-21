
for k=1:N

    disp(100*k/N);
    
    load(['Measures/' num2str(k) '.mat'])

    exc=[Measurements.Eccentricity];
    [v idx] = sort(exc,'descend');

    s=60;
    w=2*s+1;
        
%    addpath /Users/bieler/Desktop/matlab/general_functions/regionbased_seg
   
    a = imread(['zStackedYFP/' num2str(k) '.png']);
   
    N1 = size(a,1); N2 = size(a,2);    
            
    thresh = imread(['zStackedThresh/' num2str(k) '.png']);
    
    length(idx)
                        
    for p=1:length(idx)/2
        
        margin=15;
                
        %if( Measurements(idx(p)).Eccentricity > 0.8 && Measurements(idx(p)).MeanIntensity > 70)
    
        i=round(Measurements(idx(p)).Centroid(1));
        j=round(Measurements(idx(p)).Centroid(2));
    
        seli = (i-s):(i+s); %seli=mod(seli,N1);
        selj = (j-s):(j+s); %selj=mod(selj,N1);
        
        seli = abs(seli); seli(seli==0)=1;
        selj = abs(selj); selj(selj==0)=1; %reflection
        
        seli = seli -2*mod(seli.*(seli>N2),N2);
        selj = selj -2*mod(selj.*(selj>N1),N1);


        if(i>margin && j>margin && i <= N2-margin && j <= N1-margin )
         
            sub_a = a(selj,seli);
            m = thresh;
            
            px = Measurements; %erase other objects
            px=px(idx(p)).PixelIdxList; 
            m=double(m);
            m(px)=2;
            m = m==2;
                        
            m = m(selj,seli);   
            mOri=m;
            
            %is there more than 2 objects?
            labeledImage = bwlabel(m, 8);     % Label each blob so we can make measurements of it
         
            if( length(unique(labeledImage)) == 2 )
            
                        
            [y x] = find(m==1);
            
            center=zeros(1,2);
            
            center(1) = mean(x);
            center(2) = mean(y);
                
            m = bwperim(m);
                    
            [y x] = find(m==1);
            
            
            [thO rO]= cart2pol(x-center(1),y-center(2)); 
            [v ind] = sort(thO); thO=thO(ind); rO=rO(ind);
            
            chull = convhull(x,y);
            x = x(chull);   y = y(chull);
            
            x = interp1(linspace(0,1,length(x)),x,linspace(0,1,length(thO)));
            y = interp1(linspace(0,1,length(y)),y,linspace(0,1,length(thO)));
            
            if doDraw
            clf
            imagesc(m)
            hold on
            plot(x,y,'r')
            plot(center(1),center(2),'bs')
            hold off
            pause(0.5)
            end
            
            %go into polar coordinates                                    
            [th r]= cart2pol(x-center(1),y-center(2));             
            [v ind] = sort(th); th=th(ind); r=r(ind);
            
            [th,ia] = unique(th); r = r(ia);
            
            rInterp = interp1(th,r,thO);
            
            nanVal = isnan(rInterp);
            rInterp(nanVal) = 0;
                                    
            err = ( abs((rInterp - rO)) ) ;
            
            err(nanVal)=0;

            if doDraw 
                plot(thO,rInterp,'r.');
                hold on
                plot(thO,rO);
                plot(thO,err,'k')
                pause(0.5)
            end
            
            
            %look for specific pattern
            thPat = linspace(-2*pi,2*pi,200);
            pat = exp(- (thPat+0*pi/2).^2 /0.2 );% +exp(- (thPat -pi/2).^2 /0.3 );
            pat = pat.^2;
            pat = pat ./ sum(pat);
                        
            peaks = getPeakFromShape(err, pat,40,threshold,0);
            
            
%             cc = (conv(err,pat,'same'));
%             
%             [minp maxp] =peakdet(cc,0.0001)
%             
%             [v ind] = max(cc);
%             angle= thO(ind)+pi/2;
%             
%             
%             clf;
%             plot(thO,conv(err,pat,'same'))
%             hold on
%             plot(thO,err)

%%
            
            if( length(peaks) >= 2)
                
                rSep = linspace(-10,10,100);
                sepx = [rO(peaks(1)).*cos(thO(peaks(1))); rO(peaks(2)).*cos(thO(peaks(2)))] + center(1);
                sepy = [rO(peaks(1)).*sin(thO(peaks(1))); rO(peaks(2)).*sin(thO(peaks(2)))] + center(2);
                
                if doDraw
                    clf;
                    imagesc(m)
                    hold on
                    plot(sepx,sepy,'r')
                    plot(x,y,'r')
                    pause(0.6)
                end
                
                
                line = drawline([sepy(1) sepx(1)],[sepy(2) sepx(2)],size(mOri));
                sep = zeros(size(mOri));
                sep(line) = 1;
                
                 se = strel('disk',1);
                 sep = imdilate(sep,se);
                 
                 mOri = mOri .* (sep <1);  
                
                
%                 %make the separation a bit bigger
%                sepx = [(rO(peaks(1))+1).*cos(thO(peaks(1))+0.1); (rO(peaks(2))+1).*cos(thO(peaks(2))+0.1)] + center(1);
%                 sepy = [(rO(peaks(1))+1).*sin(thO(peaks(1))-0.1); (rO(peaks(2))+1).*sin(thO(peaks(2))-0.1)] + center(2);
%            
%                 line = drawline([sepy(1) sepx(1)],[sepy(2) sepx(2)],size(mOri));                                
%                 mOri(line) = 0;
%                 
%                 %make the separation a bit bigger
%                sepx = [(rO(peaks(1))+1).*cos(thO(peaks(1))-0.1); (rO(peaks(2))+1).*cos(thO(peaks(2))+0.1)] + center(1);
%                 sepy = [(rO(peaks(1))+1).*sin(thO(peaks(1))+0.1); (rO(peaks(2))+1).*sin(thO(peaks(2))-0.1)] + center(2);
%           
%                 line = drawline([sepy(1) sepx(1)],[sepy(2) sepx(2)],size(mOri));                                
%                 mOri(line) = 0;
%                 
%                 %make the separation a bit bigger
%                 sepx = [(rO(peaks(1))+1).*cos(thO(peaks(1))+0.1); (rO(peaks(2))+1).*cos(thO(peaks(2))+0.1)] + center(1);
%                 sepy = [(rO(peaks(1))+1).*sin(thO(peaks(1))+0.1); (rO(peaks(2))+1).*sin(thO(peaks(2))-0.1)] + center(2);
%                 
%                 line = drawline([sepy(1) sepx(1)],[sepy(2) sepx(2)],size(mOri));                                
%                 mOri(line) = 0;
%                 
%                 %make the separation a bit bigger
%                 sepx = [(rO(peaks(1))+1).*cos(thO(peaks(1))-0.1); (rO(peaks(2))+1).*cos(thO(peaks(2))+0.1)] + center(1);
%                 sepy = [(rO(peaks(1))+1).*sin(thO(peaks(1))-0.1); (rO(peaks(2))+1).*sin(thO(peaks(2))-0.1)] + center(2);                
% 
%                 line = drawline([sepy(1) sepx(1)],[sepy(2) sepx(2)],size(mOri));                                
%                 mOri(line) = 0;
                
               se = strel('disk',4);
               mOri = imopen(mOri,se);
                
                if doDraw
                    imagesc(mOri)
                    pause(0.5);
                end
                
                %thresh = double(thresh);
                thresh( Measurements(idx(p)).PixelIdxList ) = 0;
                thresh(selj,seli) = thresh(selj,seli) + 1*mOri;
                
            end
            
            
            
            %%
            
%             D = -bwdist(~m); 
%             L=watershed(D,8);
%             m(L == 0) = 0;
%             m = imerode(m,strel('disk',2));            
%             %imagesc(m)
%             %pause(0.1);
%             
            %seg = region_seg(sub_a, m, 80,0.2,true); %-- Run segmentation
                        
           % thresh( Measurements(idx(p)).PixelIdxList ) = 0;
           % thresh(selj,seli) = thresh(selj,seli) + seg;
            
            end
                
        end
        
    end
    
    if(doDraw)
    imagesc(thresh)
    pause(0.1);
    end
    
 
    imwrite(thresh,[saveFolder num2str(k) '.png']);

end


