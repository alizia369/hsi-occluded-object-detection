normat=mat(:,:,3:4:end-4);
[edgesum,edgecube]=getedgecube(normat);
strongedges=edgesum>1;
ns=10;
 patches=getpatchesNew(ns,strongedges); %if ns is 5 then neighbourhood is 11x11

%%% Sparse DepthMap

[ sparseMap,nneidis,meanvalpat ] = sparseDepthKLDivergence( normat,patches );

OneDresult=normaliseMap(sparseMap);
figure,imagesc(sparseMap);
%%
imgSave=true;

if(imgSave)
     imwrite(OneDresult,'.\data\snapscan b150u\Results\Depth\cube3_36bands_Depth.png');
     imwrite( ind2rgb(im2uint8(OneDresult), parula(256)), '.\data\snapscan b150u\Results\Depth\cube3_36bands_Depth-imagesc.png');
end

%%% Test getpatches

% nedges=int8(strongedges);
% 
% psz=size(nedges);
% [pr,pc]=find(nedges==1);
% edgespts=[pr,pc];
% lbl=3;
% i=1;
% totaledges=numel(pr);
% 
% %lbl=3:totaledges;
% 
% smallns=ns;
% for i=1:totaledges
%     
%    %[pr,pc]=find(nedges==1);
%   %i
%    if(numel(pr)==0)
%        break;
%    end
%   
%    fpr=pr(i);
%    fpc=pc(i);
%        
%    
%    
%    
%    
%    rmin2=fpr-smallns;%%small
%    if(rmin2<=0)
%        rmin2=1;
%        nedges( sub2ind(psz, fpr, fpc))=0;
%        continue;
%    end
%    
% 
%     
%     rmax2=fpr+smallns;
%     if(rmax2>psz(1))
%        rmax2=psz(1);
%         nedges( sub2ind(psz, fpr, fpc))=0;
%        continue;
%     end
%     
% 
%     cmin2=fpc-smallns;
%     if(cmin2<=0)
%        cmin2=1;
%         nedges( sub2ind(psz, fpr, fpc))=0;
%        continue;
%     end        
%     
%     
%     cmax2=fpc+smallns;
%     if(cmax2>=psz(2))
%        cmax2=psz(2); 
%        nedges( sub2ind(psz, fpr, fpc))=0;
%        continue;
%     end
%     
%    locnei=nedges(rmin2:rmax2,cmin2:cmax2);
%    
%    if(isempty(locnei(locnei>2)))
%     nedges(rmin2:rmax2,cmin2:cmax2)=lbl;%neighbourhood pixels
%     lbl
%     lbl=lbl+1;
%    %else
%        %nedges( sub2ind(psz, fpr, fpc))=2;%%%1=unexplored edge, 2=temp value and 3...n labels
%        %continue;
%    end     
% end
% 
% patches=nedges-2;
%  
% patches(patches<0)=0;
% 
% figure,imagesc(patches);


% %[D,I] = pdist2(edgespts,edgespts,'euclidean','Smallest',ns+1);
% [D,I] = pdist2(edgespts,edgespts,'euclidean','Smallest',totaledges);
% 
% noList=[];
% patchInd=zeros([1,totaledges]);
% j=1;
% %%
% for i=1:totaledges
%    if(~ismember(i,noList))
%        i
%         ed=D(:,i);
%         %ed
%         patchInd(j)=i;
%         j=j+1;
%         
%         pts=ed(ed<=ns+1);
%         lenpts=length(pts);
%         ptsInd=I(2:lenpts);%skipping first as it is itself start from 2
%         %ptsInd
%         noList=unique([noList,ptsInd]);
%    end
% end
% 
% 
