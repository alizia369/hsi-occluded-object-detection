
I=nendmembers;
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

        win_idx = M_idx(y-r:y+r,x-r:x+r);
        win_idx=win_idx(:);
        
        winI=I(y-r:y+r,x-r:x+r,:);
        

        winI=reshape(winI,wr,c);

        
        
        

            matCov=cov(winI);
            matCovSum=sum(matCov);
          
%%%%%SID
         
         sumWinI=sum(winI,2);
         sumWinI(sumWinI==0)=0.000001;%To avoid NAN
%          if(sumWinI==0)%This is a row %To avoid NAN
%             pw=winI;
%          else
             pw=winI./sumWinI;
%          end
         
         
         %%%%try without loop
         pwrep=repmat(pw,wr,1);
         pwirep=kron(pw,ones(wr,1));%repeat each row n times
         
%          pwrep=pwrep*matCov;
%          pwirep=pwirep*matCov;
         
         xFTerm=pwrep;
         xLogTerm=pwrep-pwirep;%pwrep./pwirep;
         
         yFTerm=pwirep;
         yLogTerm=pwirep-pwrep;%pwirep./pwrep;
         
         
         %%%%%%Nan Test
%           if(sum(isnan(xLogTerm(:)))>0)
%               disp('xLogTerm');
%             sum(isnan(xLogTerm(:)))
%           end
%           if(sum(isnan(yLogTerm(:)))>0)
%             disp('yLogTerm');
%             sum(isnan(yLogTerm(:)))
%          end
%          if(sum(isnan(matCovSum(:)))>0)
%             disp('matCovSum');
%             sum(isnan(matCovSum(:)))
%          end

        
%          smX=(sum((xFTerm.*(xLogTerm+(matCovSum+0.00005))),2))/wr ;
%          smY=(sum((yFTerm.*(yLogTerm+(matCovSum+0.00005))),2))/wr ;

         smX=(sum((xFTerm.*(xLogTerm +sum(matCovSum)+0.00005)),2))/wr;
         smY=(sum((yFTerm.*(yLogTerm +sum(matCovSum)+0.00005)),2))/wr;

         
         SID=smX+smY;
         
         tvals=SID;
            
         %%%%Testing isNaN
%          if(sum(isnan(tvals(:)))>0)
%              disp('tvals NAN');
%              sum(isnan(tvals(:)))
%          end
%          if(sum(isinf(tvals(:)))>0)
%              disp('tvals INF');
%              sum(isinf(tvals(:)))
%          end
%          

        
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

% %% %%%%%%%%%Testing
% testd=diag(Lsid);
% testimg=reshape(testd,[h w]);
% testimg=full(testimg);
% figure,imshow(testimg,[]);
% %% %%%%%%%%%%%%%%%%%Testing2
% testd=diag(Lsid,-1);
% testimg2=reshape([testd ;0],[h w]);
% testimg2=full(testimg2);
% figure,imshow(testimg2,[]);

%%% Sparse to Dense Map

%edgeMap=OneDresult;%dilatested;
std=1;
lambda=0.001;
maxBlur=3;%max(full(diag(Lsid)))+0.1;% previously 3
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

% maindiag=diag(Lsid);
% diagimg=full(reshape(maindiag,[h w]));

sparseDMap=(OneDresult.*dedgemap);
%sparseDMap=OneDresult.*dedgemap;

%sparseDMap=zeros(size(dedgemap));
% [irow,icol]=find(dedgemap~=0);%get all image indices
% imgind=[icol , irow];
% pixtopatdis = pdist2(imgind,patchind');
% [sortval,sortind]=sort(pixtopatdis,2,'ascend');
% nsortval=sortval(:,1:npat);
% nsortind=sortind(:,1:npat);
%      
% nvals=valpatchs(nsortind);
%      
% imgvals=(sum(nvals,2))./(sum(nsortval,2));
% ilinearInd = sub2ind(size(dedgemap), irow, icol);
% sparseDMap(ilinearInd)=imgvals;%reshape(imgvals,size(strongedges));
% %norresult=normaliseMap(sparseDMap);
% %figure,imshow(norresult);     
% %sparseDMap=normaliseMap(sparseDMap.^2);

% figure,imshow(sparseDMap,[]);



%%% Filter the spase defocus map to remove noise
[h,w,c]=size(normat);
sparseDMap(sparseDMap>maxBlur)=maxBlur;
sd=5; % spatial sigma
sr=0.1*maxBlur; % range sigma
[sparseDMap ~]=sparseBilateralFilter(sparseDMap,normat,sd,sr);

%figure,imshow(sparseDMap,[]);


%%% Propagate the spase defocus map based on color similarities to get full defocus map
sizeI=w*h;
constsMap=sparseDMap>0.0001;
%L=getLaplacian2(HsI,1); % calculated above
Lap=Lsid;%Lap2;%L;
D=spdiags(constsMap(:),0,sizeI,sizeI);

%%
% %%%%try pcg if you get out of memory error using the \ solver
% %%x=pcg(L+lambda*D,lambda*D*sparseDMap(:),1e-6,1000); 

%fullDMap=[];%to be removed
x=(Lap+lambda*D)\(lambda*D*sparseDMap(:));  % Solve systems of linear equations Ax = B for x or x = A\B
fullDMap=reshape(x,h,w);

figure,imshow(fullDMap,[]);





