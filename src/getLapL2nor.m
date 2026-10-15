function [ Lap ] = getLapL2nor( I,normat )
%GETLAPL2NOR Summary of this function goes here
%   Detailed explanation goes here
%I=nendmembers;
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

%vecpro=zeros([wr*wr c]);
%vecpro2=zeros([wr*wr nc]);
%%%%%%%%%%%% norms %%%%%%%%%%%
%normatrs=reshape(I,[h*w c]);
%matnorms=sqrt(sum(normatrs.^2,2));


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
          
         % sumWinI2=sum(winI2,2);
         % sumWinI2(sumWinI2==0)=0.000001;%To avoid NAN
          
           %matCov=cov(winI);
           %matCovSum=sum(matCov);
           
           
%           if (sumWinI>0)
%           pw=(winI-min(winI(:)))./(max(winI(:))-min(winI(:))); % normalise with sum
%           else
%           pw=winI;%not normalised
%           end
        
         pw=winI./sumWinI;
         pwrep=repmat(pw,wr,1);
         pwirep=kron(pw,ones(wr,1));%repeat each row n times
         
%          pw2=winI2./sumWinI2;
%          pwrep2=repmat(pw2,wr,1);
%          pwirep2=kron(pw2,ones(wr,1));%repeat each row n times
         
         
         L2dis=sqrt(sum((pwrep-pwirep).^2,2));%pwrep./pwirep;
         
         %L2dis2=sqrt(sum((pwrep2-pwirep2).^2,2));
         
         tvals=L2dis;%+L2dis2+sum(matCovSum);
        
        
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



vals=vals(1:len);
%%%%%%%%%%%%%%%%%%%%Normalise values %%%%%%%%
%vals=normaliseMap(vals);

row_idx=row_idx(1:len);
col_idx=col_idx(1:len);
Lap=sparse(row_idx,col_idx,vals,h*w,h*w);
sumL=sum(Lap,2);
Lap=spdiags(sumL(:),0,h*w,h*w)-Lap;

if(~issymmetric(Lap))% if not symmetric due to round off error do the following
Lap = (Lap+Lap.')/2; % make sure matrix is semetric without any roundoff error
end

%% Testing

% testd=diag(Lsam);
% testimg=reshape(testd,[h w]);
% testimg=full(testimg);
% figure,imshow(testimg,[]);



end

