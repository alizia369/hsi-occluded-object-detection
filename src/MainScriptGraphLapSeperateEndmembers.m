I=nendmembers;

sparseDMaps=zeros([h, w ,r]);
for i=1:r

    [ sparseDMap ] =  getProcessedSparseMap2( OneDresult,edgeMaps(:,:,i),normat );
    figure,imshow(sparseDMap);
    sparseDMaps(:,:,i)=sparseDMap;
end

%[ sparseDMap ] =  getProcessedSparseMap3( OneDresult,normat );

[h w c] = size(I);
[nh nw nc]=size(normat);
%figure,imshow(sparseDMap);
%%
%insp=find(sparseDMap~=0);
orCube=I;%normat(:,:,12:27);%normat, I
Scrible=orCube;%normat
[hs,ws,cs]=size(Scrible);
for i=1:cs %normat size can try end members as well
    %tslice=normat(:,:,i);
    insp=find(sparseDMaps(:,:,i)~=0);
    tslice=orCube(:,:,i);
    tslice(insp)=sparseDMap(insp);
    Scrible(:,:,i)=normaliseMap(tslice);
    
end

%consts_map=sum(abs(orCube-Scrible),3)>0.0001;
consts_map_all=abs(orCube-Scrible)>0.0001;
%graySc=normaliseMap(sum(Scrible,3));
graySc=Scrible;
consts_vals_all=graySc.*consts_map_all;

%%
for i=1:r

    Lap{i}=getLapL2nor( I(:,:,i),normat );
end
%%
for i=1:r
    testd=diag(Lap{i});
    testimg=reshape(testd,[h w]);
    testimg=full(testimg);
    figure,imshow(testimg,[]);
end

%%
sizeI=w*h;

%Lap=Lsam;%Lap2;%L;

lambda=0.001;
% try pcg if you get out of memory error using the \ solver+
%x=pcg(Lap+lambda*D,lambda*D*sparseDMap(:),1e-6,1000); 
fullDMaps=zeros([h,w,c]);
for i=1:r
    
   consts_map=consts_map_all(:,:,i);
   consts_vals=consts_vals_all(:,:,i);
   D=spdiags(consts_map(:),0,sizeI,sizeI);
    
x=(Lap{i}+lambda*D)\(lambda*consts_map(:).*consts_vals(:));  % Solve systems of linear equations Ax = B for x or x = A\B

fullDMaps(:,:,i)=reshape(x,h,w);

figure,imshow(fullDMaps(:,:,i),[]);

alpha=max(min(reshape(x,h,w),1),0);
figure,imshow(alpha,[]);

end
%%
Bin_imgs=zeros([h,w,c]);%r=c
for i=1:r
    Bin_imgs(:,:,i) = (imbinarize(fullDMap(:,:,i))+(i.*r));
    figure,imagesc(Bin_imgs(:,:,i));
    
end
% bw=ones([h,w]);
% for i=1:r
%     
%     bw=bw.* Bin_imgs(:,:,i);
%     
% end
% 
figure,imagesc(sum(Bin_imgs,3));