clear;
mainFolder='.\data\Data_20-7-2019\';
D = dir(mainFolder); % A is a struct ... first elements are '.' and '..' used for navigation.
 
for k = 3:length(D) % print folders and k values
 %k   
 %D(k).name  %Get the current subdirectory name
 disp([num2str(k) '. ' D(k).name]);
end

k=7; % index start from 3 monkey=4, cubesoverlap=7

currD = D(k).name;
codedir=pwd;
cd([mainFolder currD]) % change the directory (then cd('..') to get back)
worksp=[cd '\workspace.mat'];
%S=load(worksp);
load(worksp);
cd(codedir);
figure,imshow(normat(:,:,9),[]);

%%
% ns=5;
% patches=getpatches(ns,strongedges); %if ns is 5 then neighbourhood is 11x11


%%
totalp=max(patches(:));
[h w c] = size(normat);
spmap=zeros([h w]);
weights=ones([1 c]);

for i=1:totalp
   
    [rp,cp]=find(patches==i);
    %lind=find(patches==i);
    %testimg(min(rp):max(rp),min(cp):max(cp))=0;
    %testimg2(lind)=0;
    rmin=min(rp);
    rmax=max(rp);
    cmin=min(cp);
    cmax=max(cp);
    win=normat(rmin:rmax,cmin:cmax,:);
    win(:,:,1:15)=win(:,:,1:15)*1;
    win(:,:,16:24)=win(:,:,16:24)*0.5;
    win(:,:,25:34)=win(:,:,25:end)*0.2;
    
     szn=size(win);
     winI=reshape(win,szn(1)*szn(2),szn(3));
     
     covariance=cov(winI);
     dvec=diag(covariance);
     sdevc=sum(dvec);
        
     [nwr, nb]=size(winI);
     pwrep=repmat(winI,nwr,1);
     pwirep=kron(winI,ones(nwr,1));%repeat each row n times
     
     alldiff=(((pwirep.*(log(pwrep./pwirep)))) + (pwrep.*(log(pwirep./pwrep))));


        
     alldiff(isinf(alldiff))=0;
     alldiff(isnan(alldiff))=0;
        %alldiff=(subcube-npix).^2;
        
     sumvec=sum(alldiff);
        %sumvec=sum(alldiff);
        %size(sumvec)
     weights= [weights(1:6) weights(7:12) weights(13:19).*0.5 weights(20:26).*-0.5 weights(27:28).*-0.2 weights(29:31).*-0.2 weights(32:end).*-0.1];
     spmap(rmin:rmax,cmin:cmax)=sum((sumvec+(weights.*dvec')));
     %Diffcube(r(i),c(i),:)=sumvec;

    
end


figure,imagesc(spmap);
figure,imagesc(abs(spmap));
