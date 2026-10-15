function [ sparseMap,nneidis,meanvalpat ] = sparseDepthKLDivergence( normat,patches )
%SPARSEDEPTHKLDIVERGENCE Summary of this function goes here
%   Detailed explanation goes here

totalp=max(patches(:));
[h w c] = size(normat);
spmap=zeros([h w]);
weights=ones([1 c]);

nneidis=zeros(2,totalp);
meanvalpat=zeros(1,totalp);

for i=1:totalp
   
    [rp,cp]=find(patches==i);
    
     px=cp;
     py=rp;
     mpx=round(mean(px));
     mpy=round(mean(py));
     nneidis(:,i)=[mpx;mpy];
    
    
    %lind=find(patches==i);
    %testimg(min(rp):max(rp),min(cp):max(cp))=0;
    %testimg2(lind)=0;
    rmin=min(rp);
    rmax=max(rp);
    cmin=min(cp);
    cmax=max(cp);
    win=normat(rmin:rmax,cmin:cmax,:);
%     win(:,:,1:15)=win(:,:,1:15)*1;
%     win(:,:,16:24)=win(:,:,16:24)*0.5;
%     win(:,:,25:34)=win(:,:,25:34)*0.2;
    
     szn=size(win);
     winI=reshape(win,szn(1)*szn(2),szn(3));
     
     covariance=cov(winI);
     dvec=diag(covariance);
     sdevc=sum(dvec);
        
     [nwr, nb]=size(winI);
     %nwr
     %pwrep=repmat(winI,nwr,1);
     %pwirep=kron(winI,ones(nwr,1));%repeat each row n times
     
     %%%%Test only with center pixel
     tavec=squeeze(normat(mpy,mpx,:));
     avec=tavec';
     pwirep=repmat(avec,nwr,1);%rowwise
     pwrep=winI;
     %%%
     %alldiff=(((pwirep.*(log(pwrep./pwirep)))) + (pwrep.*(log(pwirep./pwrep))));
     %alldiff=(((pwirep.*(log(pwirep./pwrep)))) .* (pwrep.*(log(pwrep./pwirep))));
     alldiff=(((pwirep.*(log(pwirep./pwrep)))) + (pwrep.*(log(pwrep./pwirep))));

        
     alldiff(isinf(alldiff))=0;
     alldiff(isnan(alldiff))=0;
        %alldiff=(subcube-npix).^2;
        
     sumvec=sum(alldiff);
        %sumvec=sum(alldiff);
        %size(sumvec)
        
%     weights= [weights(1:6).*1 weights(7:12).*1 weights(13:19).*0.1 weights(20:26).*0.75 weights(27:28).*0.5 weights(29:31).*0.2 weights(32:end).*0.1];
      meanWI=mean(winI);
      Xtest = mvnrnd(meanWI,covariance,c);
      iniWeights=mvnpdf(Xtest);
      weights= (abs(log(iniWeights))+0.001)';
%      weights
%      size(meanWI)
%      size(covariance)
%      size(weights)
%      size(sumvec)
     %spmap(rmin:rmax,cmin:cmax)=sum((sumvec+(weights.*dvec')));
     %spmap(rmin:rmax,cmin:cmax)=sum((sumvec.*(weights)+dvec'.*(weights)));
     spmap(rmin:rmax,cmin:cmax)=sum((sumvec.*(weights)));
     %spmap(rmin:rmax,cmin:cmax)=sum(diff(sumvec));

     tempVal= spmap(rmin:rmax,cmin:cmax);
     meanvalpat(:,i)=mean(abs(tempVal(:)));

     %Diffcube(r(i),c(i),:)=sumvec;

    
end


% figure,imagesc(spmap);
% figure,imagesc(1-normaliseMap(spmap));

sparseMap=abs(spmap);

end

