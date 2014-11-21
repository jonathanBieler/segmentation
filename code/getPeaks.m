function maxp = getPeaks(tmp,expe,method,doDraw)
%function maxp = getPeaks(tmp,expe,method,doDraw)

NToTrack = length(tmp);

if nargin == 3
   doDraw = 0; 
end


switch method

    case 'conv'
    
    tmp(tmp==0) = mean(tmp(tmp>0));

    HP = ones(size(tmp));
    HP(1) = 0.4; HP(end) = 0.4;
    HP(2) = 0.9; HP(end-1) = 0.9;


    tmp = fft(tmp);
    tmp = tmp.*HP;

    tmp = imnorm(smooth(real(ifft(tmp)),3))';

    t = expe.t(1:NToTrack);

    f = exp(-((t-t(round(length(t)/2)))/3 ).^2)-0.02;

    pe = (conv(f,[0*tmp+tmp(1) tmp 0.8*fliplr(tmp)],'same'));

    pe = pe./(2+abs(gradient(pe)));

    pe = imnorm(pe);


    [maxp minp] = peakdet(pe,0.05);
    
    case 'diff'
        
    %%
    
%     t = linspace(0,72,140);
%     tmp = exp(-0.01*t).*( 1+cos(2*pi*t/24+6) ) +  0.4*randn(size(t)) ;
%     
%     
    idx = tmp == 0;
    if(sum(idx) > 0)
        tmp(idx) = mean(tmp(~idx)) + randn(1,sum(idx));
    end

    HG = fspecial('gaussian',100,5);
    f = imfilter(tmp,HG,'replicate');
    
    if doDraw 
       
        clfh;
        plot(tmp)
        plot(f,'r')
        pause(0.5);
        
    end
    
    f = sign( diff([f(1)-f(2) f]) );
    f = diff([f(1) f]);
    f(f>0)=0;
    f = -f;
    
    %remove peaks in zero region:
    f = f .* ~idx;
    
    maxp = find(f)';
   
    if(isempty(maxp))
        maxp = [];
    end
        
case 'fastDiff'
        
    %%
    
%     t = linspace(0,72,140);
%     tmp = exp(-0.01*t).*( 1+cos(2*pi*t/24+6) ) +  0.4*randn(size(t)) ;
%     
%     
    idx = tmp == 0;
    if(sum(idx) > 0)
        tmp(idx) = mean(tmp(~idx)) + randn(1,sum(idx));
    end

    %HG = fspecial('gaussian',100,5);
    %f = imfilter(tmp,HG,'replicate');
    f=smooth(tmp,10)';
    
    if doDraw 
       
        clfh;
        plot(tmp)
        plot(f,'r')
        pause(0.5);
        
    end
    
    f = sign( diff([f(1)-f(2) f]) );
    f = diff([f(1) f]);
    f(f>0)=0;
    f = -f;
    
    %remove peaks in zero region:
    f = f .* ~idx;
    
    maxp = find(f)';
   
    if(isempty(maxp))
        maxp = [];
    end
        
    %%
end
