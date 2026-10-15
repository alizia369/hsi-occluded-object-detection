function [ Lap ] = graphLapL2( I )
%GRAPHLAPL2 Summary of this function goes here
%   Detailed explanation goes here
eps=0.0000001;
r=1;
 [h w c] = size(I);
 [nh nw nc]=size(I);

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


%normatrs2=reshape(normat,[nh*nw nc]);
%matnorms2=sqrt(sum(normatrs2.^2,2));

tic
for y=1+r:h-r
    for x=1+r:w-r
        
        win_idx = M_idx(y-r:y+r,x-r:x+r);
        win_idx=win_idx(:);
        
        winI=I(y-r:y+r,x-r:x+r,:);
        
        %winI2=normat(y-r:y+r,x-r:x+r,:);
        
        winI=reshape(winI,wr,c);
        %winI2=reshape(winI2,wr,nc);
     %%%%%%%%% Calculating L2 for Similirity %%%%%%%%%%%%%%%%%%%%%%%%   
        
          sumWinI=sum(winI,2);
          sumWinI(sumWinI==0)=0.000001;%To avoid NAN
          if (sumWinI>0)
          pw=(winI-min(winI(:)))./(max(winI(:))-min(winI(:))); % normalise with sum
          else
          pw=winI;%not normalised
          end
        
        
         pwrep=repmat(pw,wr,1);
         pwirep=kron(pw,ones(wr,1));%repeat each row n times
         

         
         
         L2dis=sqrt(sum((pwrep-pwirep).^2,2));%pwrep./pwirep;
         
         tvals=L2dis;
        
        
        % repmat is too slow, replace repmat with indexing
        % meanI(ones(wr,1),:) <=> repmat(meanI,wr,1) 
        %winI=winI-meanI(ones(wr,1),:);
        
        %tvals=(1+winI*inv(Sigma)*winI')/wr;
        
        % win_idx(:,ones(wr,1))<=>repmat(win_idx,1,wr)
        row_idx(1+len:wr^2+len)=reshape(win_idx(:,ones(wr,1)),wr^2,1);
        
        % t(ones(wr,1),:)<=>repmat(win_idx',wr,1)
        t = win_idx';
        col_idx(1+len:wr^2+len)=reshape(t(ones(wr,1),:),wr^2,1);
        
        vals(1+len:wr^2+len)=tvals(:);
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
LL2=sparse(row_idx,col_idx,vals,h*w,h*w);
sumL=sum(LL2,2);
LL2=spdiags(sumL(:),0,h*w,h*w)-LL2;

if(~issymmetric(LL2))% if not symmetric due to round off error do the following
LL2 = (LL2+LL2.')/2; % make sure matrix is semetric without any roundoff error
end
Lap=LL2;

%% Testing

testd=diag(Lap);
testimg=reshape(testd,[h w]);
testimg=full(testimg);
figure,imshow(testimg,[]);

end

