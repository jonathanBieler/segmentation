function tmp = fillTrace(tmp)

nZ = find(tmp>0);

for i=1:length(tmp)
   if(tmp(i)==0)
      
      [v idx] = min( abs(i - nZ));
       
      tmp(i)= tmp(nZ(idx));
   end
end