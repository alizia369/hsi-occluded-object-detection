
I=nendmembers;
%[ sparseDMap ] = getProcessedSparseMap( meanvalpat,nneidis,normat );
sparseDMaps=zeros([h, w ,r]);
for i=1:r

    [ sparseDMap ] =  getProcessedSparseMap2( OneDresult,edgeMaps2(:,:,i),normat );
    figure,imshow(sparseDMap);
    sparseDMaps(:,:,i)=sparseDMap;
end

[h w c] = size(I);
[nh nw nc]=size(normat);

%%
%insp=find(sparseDMap~=0);
sparseDMap=sum(sparseDMaps,3);
insp=find(sparseDMap~=0);
orCube=I;%normat(:,:,12:27);%normat, I
Scrible=orCube;%normat
[hs,ws,cs]=size(Scrible);
for i=1:cs %normat size can try end members as well
    %tslice=normat(:,:,i);
    tslice=orCube(:,:,i);
    tslice(insp)=sparseDMap(insp);
    Scrible(:,:,i)=normaliseMap(tslice);
    
end

figure,imshow(Scrible(:,:,1));

% I=double(OriginalColorImage)/255;
% mI=double(Scrible)/255;
%consts_map=sum(abs(normat-Scrible),3)>0.001;
consts_map=sum(abs(orCube-Scrible),3)>0.001;

% if (size(I,3)==3)
graySc=normaliseMap(sum(Scrible,3));
%consts_vals=rgb2gray(mI).*consts_map;
consts_vals=graySc.*consts_map;

%%


%Lap=getLaplacianProposed(I);
%Lap=getLaplacian1(I,consts_map);
Lap=getLapL2nor( normat );

testd=diag(Lap);
testimg=reshape(testd,[h w]);
testimg=full(testimg);
figure,imshow(testimg,[]);

%%

sizeI=w*h;

%Lap=Lsam;%Lap2;%L;
D=spdiags(consts_map(:),0,sizeI,sizeI);
lambda=0.01;
% try pcg if you get out of memory error using the \ solver+
%x=pcg(Lap+lambda*D,lambda*D*sparseDMap(:),1e-6,1000); 
x=(Lap+lambda*D)\(lambda*consts_map(:).*consts_vals(:));  % Solve systems of linear equations Ax = B for x or x = A\B

fullDMap=reshape(x,h,w);

figure,imshow(fullDMap,[]);

alpha=max(min(reshape(x,h,w),1),0);
figure,imshow(alpha,[]);
%% Second way to propogate

% allDmaps=zeros([h, w ,r]);
% 
% for i=1:r
%     sparseDMap=sparseDMaps(:,:,i);
% 
%     sizeI=w*h;
%     constsMap=sparseDMap>0.0001;
% 
%     lambda=0.1;
%     %L=getLaplacian2(HsI,1); % calculated above
%     %Lap
%     D=spdiags(constsMap(:),0,sizeI,sizeI);
% 
%     % try pcg if you get out of memory error using the \ solver+
%     %x=pcg(Lap+lambda*D,lambda*D*sparseDMap(:),1e-6,1000); 
% 
%     x=(Lap+lambda*D)\(lambda*D*sparseDMap(:));  % Solve systems of linear equations Ax = B for x or x = A\B
%     fullDMap=reshape(x,h,w);
% 
%     %figure,imshow(sparseDMap,[]);
%     figure,imshow(fullDMap,[]);
%     allDmaps(:,:,i)=fullDMap;
% end