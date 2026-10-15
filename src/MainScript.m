%%%%% Load Data %%%%
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



%% Unmixing Getting given endmembers
tic
addpath(genpath('.\lib'));
cube=normat(:,:,12:27);
[h w c]=size(cube);

% Select a subset of wavelengths (17) to reduce computational time
% (Remove this line to have the same test as in the paper)
%R = R(:,:,1:10:end); 

% Run H2NMF


algo = 1; %1=rank-two NMF
%r = [];
r = 6; % number of clusters %cubesoverlap r=6.
%[IDX, C, J, sol] = hierclust2nmf(R,r,algo);
[IDX, C, J, sol] = hierclust2nmf(cube,r,algo);


%for i = 1 : size(R,3);
%for i = 1 : size(cube,3);
    %M(i,:) = reshape(R(:,:,i),sz(1)*sz(2),1);
%    M(i,:) = reshape(cube(:,:,i),sz(1)*sz(2),1);
%end 
M=reshape(cube,[h*w , c]);

H = nnlsm_blockpivot(C,M'); 

nendmembers=reshape(H',[h, w , r]);
toc
[TH,IH]=max(H);
labelimg=reshape(IH,[h,w]);
%unmixmap=testimg;
figure,imagesc(labelimg);



%% patches sparse map
%[edgesum,edgecube]=getedgecube(normat)
%strongedges=edgesum>5;
%  ns=10;
%  patches=getpatches(ns,strongedges); %if ns is 5 then neighbourhood is 11x11

%%%% Sparse DepthMap

% [ sparseMap,nneidis,meanvalpat ] = sparseDepthKLDivergence( normat,patches );
% 
% OneDresult=normaliseMap(sparseMap);
% figure,imagesc(sparseMap);

%%

% testimg=sum(nendmembers,3);
% figure,imshow(testimg,[]);
% X=zeros([h*w r]);
% for i=1:r
%     endmem=nendmembers(:,:,i);
%     X(:,i)=endmem(:);    
% end
%  class=kmeans(X,round(r/2),'MaxIter',1000);
% 
% segmap=reshape(class,h,w);
% figure,imagesc(segmap);

% endmem=nendmembers(:,:,9);
% class=kmeans(endmem(:),2);
% 
% segmap=reshape(class,h,w);
% figure,imagesc(segmap);

%%
edgeMaps=zeros([h, w ,r]);
edgeMaps2=zeros([h, w ,r]);
edgeMaps2Dir=zeros([h, w ,r]);
for i=1:r
    endmem=nendmembers(:,:,i);
    edgeMaps(:,:,i)=edge(endmem,'canny',0.15);
    [edMag,edDir]=canny(endmem);% from lib
    
    edgeMaps2(:,:,i)=edMag;
    %edgeMaps2Dir(:,:,i)=edDir;
    %figure,imshow(edgeMap);
end
% 
% figure,imshow(sum(edgeMaps,3));
% figure,imshow(sum(edgeMaps2,3));
%figure,imshow(sum(edgeMaps2Dir,3));

%% Display


for i=1:r
   figure,imshow(nendmembers(:,:,i)); 
   %figure,imshow(edgeMaps(:,:,i));
   %figure,imshow(edgeMaps2(:,:,i));
end





