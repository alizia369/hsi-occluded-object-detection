function [ sparseDMap ] = getProcessedSparseMap3( meanvalpat,normat )
%GETPROCESSEDSPARSEMAP3 Summary of this function goes here
%   Detailed explanation goes here

bspace=5;
nbands=size(normat,3);
%firstidx=1+bspace;
mididx=round((nbands-bspace*2)/2);

midBand=normat(:,:,mididx);

dedgemap=edge(midBand,'canny');

%%
%pind=find(OneDresult>0);
%allval=OneDresult(pind);

sparseDMap=meanvalpat.*dedgemap;



%%
maxBlur=1;
sparseDMap(sparseDMap>maxBlur)=maxBlur;
sd=5; % spatial sigma
sr=0.1*maxBlur; % range sigma
[sparseDMap ~]=sparseBilateralFilter(sparseDMap,normat,sd,sr);


end

