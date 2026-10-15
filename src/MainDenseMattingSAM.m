%%%%%%%%%%% should be done only once%%%%%%%% done in savecube already
% [numRows,numCols,numSpectra] = size(normat);
% scfact = mean(reshape(sqrt(sum(normat.^2,3)),numRows*numCols,1));
% normat = normat./scfact;

%% Test nor mat Data - might not be useful
%normatblr=imfilter(normat, fspecial('disk',3), 'symmetric', 'conv');
normatblr=imfilter(normat, fspecial('disk',6), 'symmetric', 'conv');
%I=normatblr;

%%

I=normat;
eps=0.0000001;
r=1;
 [h w c] = size(I);
% N = boxfilter(ones(h, w), r);
% 
% meanI_r = boxfilter(I(:, :, 1), r) ./ N;
% meanI_g = boxfilter(I(:, :, 2), r) ./ N;
% meanI_b = boxfilter(I(:, :, 3), r) ./ N;
% 
% % variance of I in each local patch
% %           rr, rg, rb
% %   Sigma = rg, gg, gb
% %           rb, gb, bb
% varI_rr = boxfilter(I(:, :, 1).*I(:, :, 1), r) ./ N - meanI_r .*  meanI_r;
% varI_rg = boxfilter(I(:, :, 1).*I(:, :, 2), r) ./ N - meanI_r .*  meanI_g;
% varI_rb = boxfilter(I(:, :, 1).*I(:, :, 3), r) ./ N - meanI_r .*  meanI_b;
% varI_gg = boxfilter(I(:, :, 2).*I(:, :, 2), r) ./ N - meanI_g .*  meanI_g;
% varI_gb = boxfilter(I(:, :, 2).*I(:, :, 3), r) ./ N - meanI_g .*  meanI_b;
% varI_bb = boxfilter(I(:, :, 3).*I(:, :, 3), r) ./ N - meanI_b .*  meanI_b;

wr = (2*r+1) * (2*r+1);
cc = zeros(h,w);
tlen=sum(sum(1-cc(r+1:end-r,r+1:end-r)))*(wr^2);
M_idx = reshape([1:h*w], h, w);
vals=zeros(tlen,1);
vals2=zeros(tlen,1);
len = 0;

vecpro=zeros([wr*wr c]);
%%%%%%%%%%%% norms %%%%%%%%%%%
normatrs=reshape(normat,[h*w c]);
matnorms=sqrt(sum(normatrs.^2,2));

tic
for y=1+r:h-r
    for x=1+r:w-r
%         Sigma = [varI_rr(y, x), varI_rg(y, x), varI_rb(y, x);
%             varI_rg(y, x), varI_gg(y, x), varI_gb(y, x);
%             varI_rb(y, x), varI_gb(y, x), varI_bb(y, x)];
%         meanI = [meanI_r(y,x), meanI_g(y,x), meanI_b(y,x)];
     
%        Sigma = Sigma + eps * eye(3);
        win_idx = M_idx(y-r:y+r,x-r:x+r);
        win_idx=win_idx(:);
        
        winI=I(y-r:y+r,x-r:x+r,:);
        
        winI=reshape(winI,wr,c);
        
     %%%%%%%%% Calculating Spectral Angular Mapper for Similirity %%%%%%%%%%%%%%%%%%%%%%%%   
        
        for i=1:wr
            tempvals=winI.*winI(i,:);
            vecpro(((i-1)*wr)+1:wr*i,:)=tempvals;
        end
        %size(vecpro)
        sumVec=sum(vecpro,2);
        sumVecRs=reshape(sumVec(:),[wr wr]);
        %size(sumVecRs)
        
        curnorms=matnorms(win_idx);
        allcombnorms=curnorms*curnorms';
        %size(allcombnorms)
        
        divVec=sumVecRs./allcombnorms;
        
        %sad=round(acos(divVec),6); % Spectral Angular Distance for all Vectors for sad 0 is same and 1 is different.
        sad=acos(divVec);% Spectral Angular Distance for all Vectors for sad 0 is same and 1 is different.
       
        
        %%%%%% Test %%%%%
        %sad=sad+diag(sum(sad,2)); %%%put sum instead of zeros in diagnoal 
        
        similirity=1-sad; %similirity values filped 1 for same and 0 for different
        
        
        
        
        % repmat is too slow, replace repmat with indexing
        % meanI(ones(wr,1),:) <=> repmat(meanI,wr,1) 
        %winI=winI-meanI(ones(wr,1),:);
        
        %tvals=(1+winI*inv(Sigma)*winI')/wr;
        
        % win_idx(:,ones(wr,1))<=>repmat(win_idx,1,wr)
        row_idx(1+len:wr^2+len)=reshape(win_idx(:,ones(wr,1)),wr^2,1);
        
        % t(ones(wr,1),:)<=>repmat(win_idx',wr,1)
        t = win_idx';
        col_idx(1+len:wr^2+len)=reshape(t(ones(wr,1),:),wr^2,1);
        
        vals(1+len:wr^2+len)=sad(:);%tvals(:);
        vals2(1+len:wr^2+len)=similirity(:);
        
        len=len+wr^2;
    end
end

toc
%%
vals=vals(1:len);
%%%%%%%%%%%%%%%%%%%%Normalise values %%%%%%%%
%vals=normaliseMap(vals);

row_idx=row_idx(1:len);
col_idx=col_idx(1:len);
Lsam=sparse(row_idx,col_idx,vals,h*w,h*w);
sumL=sum(Lsam,2);
Lsam=spdiags(sumL(:),0,h*w,h*w)-Lsam;

if(~issymmetric(Lsam))% if not symmetric due to round off error do the following
Lsam = (Lsam+Lsam.')/2; % make sure matrix is semetric without any roundoff error
end
%% %%% values 2
vals2=vals2(1:len);
row_idx=row_idx(1:len);
col_idx=col_idx(1:len);
Lflipsam=sparse(row_idx,col_idx,vals2,h*w,h*w);
sumL2=sum(Lflipsam,2);
Lflipsam=spdiags(sumL2(:),0,h*w,h*w)-Lflipsam;

if(~issymmetric(Lflipsam))% if not symmetric due to round off error do the following
Lflipsam = (Lflipsam+Lflipsam.')/2; % make sure matrix is semetric without any roundoff error
end


%% Testing

testd=diag(Lsam);
testimg=reshape(testd,[h w]);
testimg=full(testimg);
figure,imshow(testimg,[]);

% testd2=diag(Lflipsam);
% testimg2=reshape(testd2,[h w]);
% testimg2=full(testimg2);
% figure,imshow(testimg2,[]);

%% Sparse to Dense Map

edgeMap=OneDresult;%dilatested;
std=1;
lambda=0.001;
maxBlur=1;
%sparseDMap=dilatested; %%% Can be changed
bspace=5;
nbands=size(normat,3);
%firstidx=1+bspace;
mididx=round((nbands-bspace*2)/2);

midBand=normat(:,:,mididx);

dedgemap=edge(midBand,'canny');

%%%%%%%%%%%%%%%%Morphology%%%%%%%%%%%
% se90 = strel('line', 3, 90);
% se0 = strel('line', 3, 0);
% 
% dedgemap = imdilate(dedgemap, [se90 se0]);
% figure, imshow(BWsdil), title('dilated gradient mask');


%figure,imshow(dedgemap);

%% Generate the sparse defocus map%%


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
figure,imshow(sparseDMap,[]);



%% Filter the spase defocus map to remove noise
[h,w,c]=size(normat);
sparseDMap(sparseDMap>maxBlur)=maxBlur;
sd=5; % spatial sigma
sr=0.1*maxBlur; % range sigma
[sparseDMap ~]=sparseBilateralFilter(sparseDMap,normat,sd,sr);

%% Propagate the spase defocus map based on color similarities to get full defocus map
sizeI=w*h;
constsMap=sparseDMap>0.0001;
%L=getLaplacian2(HsI,1); % calculated above
Lap=Lsam;%Lap2;%L;
D=spdiags(constsMap(:),0,sizeI,sizeI);

% try pcg if you get out of memory error using the \ solver
% x=pcg(L+lambda*D,lambda*D*sparseDMap(:),1e-6,1000); 
x=(Lap+lambda*D)\(lambda*D*sparseDMap(:));  % Solve systems of linear equations Ax = B for x or x = A\B
fullDMap=reshape(x,h,w);

figure,imshow(fullDMap,[]);

%% Propogate Based on Distance


valpatchs=meanvalpat;
patchind=nneidis;%%%% nneidis from main Script
npat=10;

result=zeros(size(strongedges));
[irow,icol]=find(result>-1);%get all image indices
imgind=[icol , irow];
pixtopatdis = pdist2(imgind,patchind');
[sortval,sortind]=sort(pixtopatdis,2,'ascend');
nsortval=sortval(:,1:npat);
nsortind=sortind(:,1:npat);
     
nvals=valpatchs(nsortind);
     
imgvals=(sum(nvals,2))./(sum(nsortval,2));

result=reshape(imgvals,size(strongedges));
norresult=normaliseMap(result);
figure,imshow(norresult);     

%%

blur1 = imgaussfilt(norresult,6);
figure,imshow(blur1);


