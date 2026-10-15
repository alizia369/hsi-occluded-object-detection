
%%

I=nendmembers;
eps=0.0000001;
r=1;
 [h w c] = size(I);


wr = (2*r+1) * (2*r+1);
cc = zeros(h,w);
tlen=sum(sum(1-cc(r+1:end-r,r+1:end-r)))*(wr^2);
M_idx = reshape([1:h*w], h, w); %#ok<NBRAK>
vals=zeros(tlen,1);
vals2=zeros(tlen,1);
len = 0;

%vecpro=zeros([wr*wr c]);%81x34

vecpro=zeros([wr*wr 1]);

%%%%Spatial Distance
valpatchs=meanvalpat;
patchind=nneidis;%%%% nneidis from main Script
npat=10;
pindin=patchind';
%%%%%
xysize=[h w];



%%%%%%%%%%Speedup
wrsq=wr*wr;

xLogTerm=zeros([wr*wr c]);%pw./pw(i,:);
xFTerm= zeros([wr*wr c]);%pw
            
yLogTerm=zeros([wr*wr c]);%pw(i,:)./pw;
yFTerm=zeros([wr*wr c]);%pw(i,:);

FullLog=zeros([2*wr*wr c]);
FullLogRes=zeros([2*wr*wr c]);


tic
for y=1+r:h-r
    for x=1+r:w-r

        win_idx = M_idx(y-r:y+r,x-r:x+r);
        win_idx=win_idx(:);
        
        winI=I(y-r:y+r,x-r:x+r,:);
        

        winI=reshape(winI,wr,c);
        

        

        
%%%%%SID
         
         sumWinI=sum(winI,2);
         pw=winI./sumWinI;
         
         
         %%%%try without loop
         pwrep=repmat(pw,wr,1);
         pwirep=kron(pw,ones(wr,1));%repeat each row n times
         
         xFTerm=pwrep;
         xLogTerm=pwrep./pwirep;
         
         yFTerm=pwirep;
         yLogTerm=pwirep./pwrep;
         
         
         %%%%%%
         
         
         FullLog(1:wrsq,:)=xLogTerm;
         FullLog(wrsq+1:end,:)=yLogTerm;
         
         FullLogRes=log(FullLog);
         
         smX=(sum(xFTerm.*FullLogRes(1:wrsq,:),2));
         smY=(sum(yFTerm.*FullLogRes(wrsq+1:end,:),2));
         SID=smX+smY;
         
         tvals=SID;
            

        
        % win_idx(:,ones(wr,1))<=>repmat(win_idx,1,wr)
        %row_idx(1+len:wr^2+len)=reshape(win_idx(:,ones(wr,1)),wr^2,1);
        row_idx(1+len:wr^2+len)=reshape(repmat(win_idx,1,wr),wr^2,1);
        
        % t(ones(wr,1),:)<=>repmat(win_idx',wr,1)
        t = win_idx';
        %col_idx(1+len:wr^2+len)=reshape(t(ones(wr,1),:),wr^2,1);
        col_idx(1+len:wr^2+len)=reshape(repmat(t,wr,1),wr^2,1);
        
        vals(1+len:wr^2+len)=tvals(:);
        len=len+wr^2;
        
         
    end
end

toc
%% %% Remove NAN and inf
vals(isnan(vals))=0;
vals(isinf(vals))=0;
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

