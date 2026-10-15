
%%% Computing GraphLapacian on endmembers from MainScript

I=nendmembers;
eps=0.0000001;
r=1;
 [h w c] = size(I);
 [nh nw nc]=size(normat);

wr = (2*r+1) * (2*r+1);
cc = zeros(h,w);
tlen=sum(sum(1-cc(r+1:end-r,r+1:end-r)))*(wr^2);
M_idx = reshape([1:h*w], h, w);
vals=zeros(tlen,1);
vals2=zeros(tlen,1);
len = 0;

vecpro=zeros([wr*wr c]);
vecpro2=zeros([wr*wr nc]);
%%%%%%%%%%%% norms %%%%%%%%%%%
normatrs=reshape(I,[h*w c]);
matnorms=sqrt(sum(normatrs.^2,2));


normatrs2=reshape(normat,[nh*nw nc]);
matnorms2=sqrt(sum(normatrs2.^2,2));

tic
for y=1+r:h-r
    for x=1+r:w-r
        
        win_idx = M_idx(y-r:y+r,x-r:x+r);
        win_idx=win_idx(:);
        
        winI=I(y-r:y+r,x-r:x+r,:);
        
        winI2=normat(y-r:y+r,x-r:x+r,:);
        
        winI=reshape(winI,wr,c);
        winI2=reshape(winI2,wr,nc);
     %%%%%%%%% Calculating Spectral Angular Mapper for Similirity %%%%%%%%%%%%%%%%%%%%%%%%   
        
        for i=1:wr
            tempvals=winI.*winI(i,:);
            vecpro(((i-1)*wr)+1:wr*i,:)=tempvals;
            
            tempvals2=winI2.*winI2(i,:);
            vecpro2(((i-1)*wr)+1:wr*i,:)=tempvals2;
        end
        %size(vecpro)
        sumVec=sum(vecpro,2);
        sumVecRs=reshape(sumVec(:),[wr wr]);
        
        sumVec2=sum(vecpro2,2);
        sumVecRs2=reshape(sumVec2(:),[wr wr]);
        %size(sumVecRs)
        
        curnorms=matnorms(win_idx);
        allcombnorms=curnorms*curnorms';
        
        curnorms2=matnorms2(win_idx);
        allcombnorms2=curnorms2*curnorms2';
        %size(allcombnorms)
        
        divVec=sumVecRs./allcombnorms;
        
        divVec2=sumVecRs2./allcombnorms2;
        
        %sad=round(acos(divVec),6); % Spectral Angular Distance for all Vectors for sad 0 is same and 1 is different.
        sad=acos(divVec);% Spectral Angular Distance for all Vectors for sad 0 is same and 1 is different.
       
        sad2=acos(divVec2);
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
        
        vals(1+len:wr^2+len)=sad(:).*sad2(:);%tvals(:);
        %vals2(1+len:wr^2+len)=similirity(:);
        
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
Lsam=sparse(row_idx,col_idx,vals,h*w,h*w);
sumL=sum(Lsam,2);
Lsam=spdiags(sumL(:),0,h*w,h*w)-Lsam;

if(~issymmetric(Lsam))% if not symmetric due to round off error do the following
Lsam = (Lsam+Lsam.')/2; % make sure matrix is semetric without any roundoff error
end
%% %%% values 2
% vals2=vals2(1:len);
% row_idx=row_idx(1:len);
% col_idx=col_idx(1:len);
% Lflipsam=sparse(row_idx,col_idx,vals2,h*w,h*w);
% sumL2=sum(Lflipsam,2);
% Lflipsam=spdiags(sumL2(:),0,h*w,h*w)-Lflipsam;
% 
% if(~issymmetric(Lflipsam))% if not symmetric due to round off error do the following
% Lflipsam = (Lflipsam+Lflipsam.')/2; % make sure matrix is semetric without any roundoff error
% end


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
figure,imshow(sparseDMap,[]);



%% Filter the spase defocus map to remove noise
[h,w,c]=size(normat);
sparseDMap(sparseDMap>maxBlur)=maxBlur;
sd=5; % spatial sigma
sr=0.1*maxBlur; % range sigma
[sparseDMap ~]=sparseBilateralFilter(sparseDMap,normat,sd,sr);

%% Propagate the spase defocus map based on color similarities to get full defocus map
sizeI=w*h;
lambda=0.01;
constsMap=sparseDMap>0.0001;
%L=getLaplacian2(HsI,1); % calculated above
Lap=Lsam;%Lap2;%L;
D=spdiags(constsMap(:),0,sizeI,sizeI);

% try pcg if you get out of memory error using the \ solver+
% x=pcg(L+lambda*D,lambda*D*sparseDMap(:),1e-6,1000); 
x=(Lap+lambda*D)\(lambda*D*sparseDMap(:));  % Solve systems of linear equations Ax = B for x or x = A\B
fullDMap=reshape(x,h,w);

figure,imshow(fullDMap,[]);