function [ patches ] = getpatchesNew(  smallns,strongedges )
%GETPATCHESNEW Summary of this function goes here
%   Detailed explanation goes here
nedges=int32(strongedges);

psz=size(nedges);
[pr,pc]=find(nedges==1);

lbl=3;

totaledges=numel(pr);

for i=1:totaledges
    
   %[pr,pc]=find(nedges==1);
  %i
   if(numel(pr)==0)
       break;
   end
  
   fpr=pr(i);
   fpc=pc(i);
       
   
   
   
   
   rmin2=fpr-smallns;%%small
   if(rmin2<=0)
       rmin2=1;
       nedges( sub2ind(psz, fpr, fpc))=0;
       continue;
   end
   

    
    rmax2=fpr+smallns;
    if(rmax2>psz(1))
       rmax2=psz(1);
        nedges( sub2ind(psz, fpr, fpc))=0;
       continue;
    end
    

    cmin2=fpc-smallns;
    if(cmin2<=0)
       cmin2=1;
        nedges( sub2ind(psz, fpr, fpc))=0;
       continue;
    end        
    
    
    cmax2=fpc+smallns;
    if(cmax2>=psz(2))
       cmax2=psz(2); 
       nedges( sub2ind(psz, fpr, fpc))=0;
       continue;
    end
    
   locnei=nedges(rmin2:rmax2,cmin2:cmax2);
   
   if(isempty(locnei(locnei>2)))
    nedges(rmin2:rmax2,cmin2:cmax2)=lbl;%neighbourhood pixels
    %lbl
    lbl=lbl+1;
   %else
       %nedges( sub2ind(psz, fpr, fpc))=2;%%%1=unexplored edge, 2=temp value and 3...n labels
       %continue;
   end     
end

patches=nedges-2;
 
patches(patches<0)=0;

%figure,imagesc(patches);


end

