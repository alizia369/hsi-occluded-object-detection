function [ fullDMap ] = graphSolveLap( Lap,I,meanvalpat,nneidis )
%GRAPHSOLVELAP Summary of this function goes here
%   Detailed explanation goes here

% single band
%% Sparse to Dense Map

std=1;
lambda=0.001;
maxBlur=1;
%sparseDMap=dilatested; %%% Can be changed
%BW = imbinarize(I);

[sh,sw,sc]=size(I);
dedgemaps=zeros(size(I));
for i=1:sc
dedgemaps(:,:,i)=edge(I(:,:,i),'canny');
end
dedgemap=sum(dedgemaps,3)>0;
figure,imshow(dedgemap);

%% Generate the sparse defocus map%%

valpatchs=meanvalpat;%%%%  from main Script
patchind=nneidis;%%%% nneidis from main Script
npat=10;

sparseDMap=zeros(size(dedgemap));

[irow,icol]=find(dedgemap~=0);%get all image indices
imgind=[icol , irow];
pixtopatdis = pdist2(imgind,patchind');
[sortval,sortind]=sort(pixtopatdis,2,'ascend');
nsortval=sortval(:,1:npat);
nsortind=sortind(:,1:npat);
     
nvals=valpatchs(nsortind);
     
imgvals=(sum(nvals,2))./(sum(nsortval,2));
ilinearInd = sub2ind(size(dedgemap), irow, icol);
sparseDMap(ilinearInd)=imgvals;%reshape(imgvals,size(strongedges));
%norresult=normaliseMap(sparseDMap);
%figure,imshow(norresult);     
%sparseDMap=normaliseMap(sparseDMap.^2);
sparseDMap=normaliseMap(sparseDMap);
%figure,imshow(sparseDMap);

%% Filter the spase defocus map to remove noise
%[h,w,c]=size(normat);
[h,w,c]=size(I);
sparseDMap(sparseDMap>maxBlur)=maxBlur;
sd=5; % spatial sigma
sr=0.1*maxBlur; % range sigma
[sparseDMap ~]=sparseBilateralFilter(sparseDMap,I,sd,sr);

%%  alpha matting parameters

insp=find(sparseDMap~=0);
Scrible=I;%normat
[hs,ws,cs]=size(Scrible);
for i=1:cs %normat size can try end members as well
    %tslice=normat(:,:,i);
    tslice=I(:,:,i);
    tslice(insp)=sparseDMap(insp);
    Scrible(:,:,i)=normaliseMap(tslice);
    
end

%figure,imshow(Scrible(:,:,1));

% I=double(OriginalColorImage)/255;
% mI=double(Scrible)/255;
%consts_map=sum(abs(normat-Scrible),3)>0.001;
consts_map=sum(abs(I-Scrible),3)>0.001;

% if (size(I,3)==3)
graySc=normaliseMap(sum(Scrible,3));
%consts_vals=rgb2gray(mI).*consts_map;
consts_vals=graySc.*consts_map;
% end
% if (size(I,3)==1)

sizeI=w*h;

%Lap=LL2;%Lap2;%L;
D=spdiags(consts_map(:),0,sizeI,sizeI);
lambda=0.09;%0.01;
% try pcg if you get out of memory error using the \ solver+
%x=pcg(Lap+lambda*D,lambda*D*sparseDMap(:),1e-6,1000); 
x=(Lap+lambda*D)\(lambda*consts_map(:).*consts_vals(:));  % Solve systems of linear equations Ax = B for x or x = A\B

%fullDMap=reshape(x,h,w);

%figure,imshow(fullDMap,[]);

alpha=max(min(reshape(x,h,w),1),0);
fullDMap=alpha;
%figure,imshow(alpha,[]);

end

