%[seli selj] = getNeiInd(i,j,s,N1,N2)
function  [seli selj] = getNeiInd(i,j,s,N1,N2)
    

    seli = (i-s):(i+s); 
    selj = (j-s):(j+s); 

    seli = abs(seli); seli(seli==0)=1;
    selj = abs(selj); selj(selj==0)=1; %reflection

    seli = seli -2*mod(seli.*(seli>N2),N2);
    selj = selj -2*mod(selj.*(selj>N1),N1);

   