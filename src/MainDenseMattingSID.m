%%%%%%%%%%% should be done only once%%%%%%%% done in savecube already
% [numRows,numCols,numSpectra] = size(normat);
% scfact = mean(reshape(sqrt(sum(normat.^2,3)),numRows*numCols,1));
% normat = normat./scfact;

%%

I=normat;
eps=0.0000001;
r=1;
 [h w c] = size(I);


wr = (2*r+1) * (2*r+1);
cc = zeros(h,w);
tlen=sum(sum(1-cc(r+1:end-r,r+1:end-r)))*(wr^2);
M_idx = reshape([1:h*w], h, w);
vals=zeros(tlen,1);
vals2=zeros(tlen,1);
len = 0;

%vecpro=zeros([wr*wr c]);%81x34
vecpro=zeros([wr*wr 1]);

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
%%%%%%%%%%%%%%%%%% Testing Diffferent Methods %%%%%%%%%%%%%%%%%%%%%%%        
        %%%%%%%%covariance %%%%%%%%%% alright
%         sigmaw=cov(winI);
%         sigmaw=sigmaw + eps * eye(c);
        %%%%%%%%% corellation coefficent %%%% not good
%         sigmaw=corrcoef(winI);
%         sigmaw=sigmaw + eps * eye(c);
       

       %%%%%%%%%%%%%%%%%%% Spectral_Information Divergence(SID)
% pX=X./sum(X);
% pY=Y./sum(Y);
% smX=sum(pX.*log(pX./pY));
% smY=sum(pY.*log(pY./pX));
% 
% SID=smX+smY;
         %winI=gradient(winI);% need to be removed just to test
           %local normalise windows
%         maxWin=max(winI(:));
%         minWin=min(winI(:));
%         winI=(winI-maxWin)./ (maxWin-minWin);
         
         
         sumWinI=sum(winI,2);
         pw=winI./sumWinI;
         
         for i=1:wr
            smX=sum(pw.*log(pw./pw(i,:)),2); %size(smX)=9x1
            smY=sum(pw(i,:).*log(pw(i,:)./pw),2);
            SID=smX+smY;%size(SID)=9x1
%             tempvals=winI.*winI(i,:);
            vecpro(((i-1)*wr)+1:wr*i)=SID;
         end
         
         tvals=vecpro;
            
%%%%%%%%%%%%% Post calculation for some methods like cov, xcov

        %%%%%%% mean
%         meanI=mean(winI);
%         winI=winI-meanI(ones(wr,1),:);
        
        %%%%%%%%%%tvalue
        
%         tvals=(1+winI*inv(sigmaw)*winI')/wr; % result= 9x9 matrix
         
        
                
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End Testing Methods         
        
        % win_idx(:,ones(wr,1))<=>repmat(win_idx,1,wr)
        row_idx(1+len:wr^2+len)=reshape(win_idx(:,ones(wr,1)),wr^2,1);
        
        % t(ones(wr,1),:)<=>repmat(win_idx',wr,1)
        t = win_idx';
        col_idx(1+len:wr^2+len)=reshape(t(ones(wr,1),:),wr^2,1);
        
        vals(1+len:wr^2+len)=tvals(:);
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
Lsid=sparse(row_idx,col_idx,vals,h*w,h*w);
sumL=sum(Lsid,2);
Lsid=spdiags(sumL(:),0,h*w,h*w)-Lsid;

if(~issymmetric(Lsid))% if not symmetric due to round off error do the following
Lsid = (Lsid+Lsid.')/2; % make sure matrix is semetric without any roundoff error
end

%% %%%%%%%%%Testing
testd=diag(Lsid);
testimg=reshape(testd,[h w]);
testimg=full(testimg);
figure,imshow(testimg,[]);
%% %%%%%%%%%%%%%%%%%Testing2
testd=diag(Lsid,-1);
testimg=reshape([testd ;0],[h w]);
testimg=full(testimg);
figure,imshow(testimg,[]);

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

patchind=nneidis;%%%% nneidis from main Script
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
Lap=Lsid;%Lap2;%L;
D=spdiags(constsMap(:),0,sizeI,sizeI);

% try pcg if you get out of memory error using the \ solver
% x=pcg(L+lambda*D,lambda*D*sparseDMap(:),1e-6,1000); 
x=(Lap+lambda*D)\(lambda*D*sparseDMap(:));  % Solve systems of linear equations Ax = B for x or x = A\B
fullDMap=reshape(x,h,w);

figure,imshow(fullDMap,[]);

