function [ Lap ] = getLaplacianProposed( normat )
%GETLAPLACIANPROPOSED Summary of this function goes here
%   Detailed explanation goes here

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
         
         


         smX=(sum((xFTerm.*(xLogTerm +sum(matCovSum)+0.00005)),2))/wr;
         smY=(sum((yFTerm.*(yLogTerm +sum(matCovSum)+0.00005)),2))/wr;

         
         SID=smX+smY;
         
         tvals=SID;
            


        
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



end

