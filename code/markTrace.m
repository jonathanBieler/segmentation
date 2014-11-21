function trace = markTrace(x)

trace = double(x > 0);

d = diff([trace trace(end)]);

ind = find(d~=0);
%ind

if(length(ind)==1 && trace(1) == 0 ) %black at begining      
    trace(ind) = -1;    
end

if(length(ind)==1 && trace(1) == 1 ) %black at the end      
    trace(ind) = -2;    
end

if(length(ind)==2 && trace(1) == 0) %one block
    
    trace(ind(1)) = -1;
    trace(ind(2)+1) = -2;
end

if(length(ind)==2 && trace(1) == 1) %one gap
    
    trace(ind(1)+1:ind(2))=-3; 
end

if(length(ind)>2 &&  trace(1) == 0 && trace(end) == 0  ) 
    
    trace(ind(1)) = -1;
    trace(ind(end)+1) = -2;
    ind = ind(2:end-1);
    
    for i=1:2:length(ind)       
        trace(ind(i)+1:ind(i+1))=-3;        
    end
    
end

if(length(ind)>2 &&  trace(1) == 1 && trace(end) == 0  ) 
        
    trace(ind(end)+1) = -2;
    ind = ind(1:end-1);
    
    for i=1:2:length(ind)       
        trace(ind(i)+1:ind(i+1))=-3;        
    end
    
end

if(length(ind)>2 &&  trace(1) == 0 && trace(end) == 1  ) 
    
    trace(ind(1)) = -1;
 
    ind = ind(2:end);
    
    for i=1:2:length(ind)       
        trace(ind(i)+1:ind(i+1))=-3;        
    end
    
end

if(length(ind)>2 &&  trace(1) == 1 && trace(end) == 1  ) 
    
    
    for i=1:2:length(ind)       
        trace(ind(i)+1:ind(i+1))=-3;        
    end
    
end

% clf; hold on
% plot(trace)
% plot(d,'r')
% pause(1);