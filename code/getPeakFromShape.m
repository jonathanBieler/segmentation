%function peaks = getPeakFromShape(y, pat,peakSize,thresh,doDraw)
function peaks = getPeakFromShape(y, pat,peakSize,thresh,doDraw)

y = y(:);
N=length(y);
x = 1:N;


c = conv([y; y; y],pat,'same')';
c = c((N+1):(2*N));

peaks = [];

while 1

    [v ind] = max(c);
    
    if doDraw
        clf;
       plot(c)
        pause(1)
    end
   
    if v > thresh

        peaks = [peaks ind];
        xx = min(abs(N-x+ind),abs(x-ind));%circular x        
        
        eraseFun = 1-exp( -xx.^2 / (peakSize)^2 );

        c = c .*eraseFun;
    else
        break;
    
    end

end


if doDraw
    clf; hold on
    plot(x,y)
    plot(x(peaks),y(peaks),'r.')
    drawnow
end

