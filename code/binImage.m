function a = binImage(a,p)

     a = double(a);

     [m,n]=size(a); %a is the original aatrix

     a=sum( reshape(a,p,[]) ,1 );
     a=reshape(a,m/p,[]).'; 

     a=sum( reshape(a,p,[]) ,1);
     a=reshape(a,n/p,[]).'; 