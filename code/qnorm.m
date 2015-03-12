function x = qnorm(x,q)


    x = x-quantile(x(:),q);   
    q = quantile(x(:),1-q);
    
    if(q>0)
        x = x/q;
    end
    
    if( sum(isnan(x(:))))
       'wtf' 
    end
