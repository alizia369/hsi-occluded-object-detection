function [ nDiffMap ] = SumofDifferenceMap( normat,edgesum,strongedge,neighbourhoodsize )
%SUMOFDIFFERENCE Summary of this function goes here
%   Detailed explanation goes here


if nargin < 3
    strongedge=5;%10;%10
    neighbourhoodsize=8;%8
    %histogrambins=15;%15

    
end

tic
sz=size(normat);

ns=neighbourhoodsize;%neighbourhood size

[r,c]=find(edgesum>=strongedge); % 5 is fine as well


Diffcube=zeros(sz(1),sz(2),sz(3));


    
   % slice = normat(:,:,1);
    
    %  tic
 
    for i=1:length(r)
        
        rmin=r(i)-ns;
        if(rmin<=0)
        rmin=1;
        end
        
        rmax=r(i)+ns;
        if(rmax>sz(1))
        rmax=sz(1);
        end
        
        cmin=c(i)-ns;
        if(cmin<=0)
        cmin=1;
        end        
        
        cmax=c(i)+ns;
        if(cmax>=sz(2))
           cmax=sz(2); 
        end
        
        
        %npix=slice(rmin:rmax,cmin:cmax);%neighbourhood pixels
        
        
        %nsz=size(npix(:));
        %allerr=[];
        %alllog=[];
        alldiff=[];
%             for j=rmin:rmax
%            
%                 for k=cmin:cmax
%                     
%                 avec=squeeze(normat(r(i),c(i),:));
%                 bvec=squeeze(normat(r(j),c(k),:));
% 
%                diffab=((avec-bvec).^2);
%                
%                alldiff=[alldiff;diffab'];
%                
%                
%                 end
%             end
        npix=normat(rmin:rmax,cmin:cmax,:);
        szn=size(npix);
%         tones=ones(szn);
%         subcube=ones(szn);
        npix2=reshape(npix,szn(1)*szn(2),szn(3));
        %avec=mean(npix2);
        tavec=squeeze(normat(r(i),c(i),:));
        avec=tavec';
        %npix2in=npix2';
        
%         for j=1:szn(3)
%            
%            subcube(:,:,j) =tones(:,:,j).*avec(j);
%             
%         end
        
       % size(npix)
       % size(subcube)
        %submat=repmat(avec,1,size(npix2in,2));
        %submat=diag(avec)*ones(size(npix2in));
        
        submat=repmat(avec,size(npix2,1),1);%rowwise
        %size(submat)
        %size(npix2)
        
        alldiff=(abs(npix2-submat));
        %alldiff=(subcube-npix).^2;
        
        sumvec=sum(alldiff);
        %size(sumvec)
        Diffcube(r(i),c(i),:)=sumvec;
        %i
    end
   
    sqrld2=Diffcube.^2;
    sumld2=sum(sqrld2,3);
    DiffMap=sqrt(sumld2);
    
   [ nDiffMap ] = normaliseMap( DiffMap );
    
    
toc
end

