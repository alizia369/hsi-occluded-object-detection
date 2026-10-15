% mainFolder='.\data\snapscan b150u\final\outside\';
mainFolder='.\data\Occluded Onject Detection final\Atof Camera\objs1\';
resultFolder=[ mainFolder , 'Results\'];
D = dir([mainFolder , '*.mat']); % A is a struct ... first elements are '.' and '..' used for navigation.
 
for i = 1:length(D) % print folders and k values
 
 disp([num2str(i) '. ' D(i).name]);
end

k=2;
name=D(k).name;
load([mainFolder , name]);


%finalMat=datacube;% need to be commented only for cave dataset
finalMat=normat;%(:,:,19:end-11);% need to be commented


mat=finalMat;

% for i=1:size(mat,3)
%    mat(:,:,i)=imadjust(mat(:,:,i)); 
% end

figure,imshow(mat(:,:,12));

%%
normat=mat;%mat(:,:,3:3:end-3);%mat(:,:,3:4:end-4);%mat;

[h,w,bands]=size(normat);
[edgesum,edgecube]=getedgecube(normat);
strongedges=edgesum>1;
ns=5;
 patches=getpatchesNew(ns,strongedges); %if ns is 5 then neighbourhood is 11x11

%%% Sparse DepthMap

[ sparseMap,nneidis,meanvalpat ] = sparseDepthKLDivergence( normat,patches );
%[ sparseMap ] = SumofDifferenceMap( normat,edgesum,strongedges,ns );%testing

OneDresult=normaliseMap(sparseMap);
figure,imagesc(sparseMap);
%% 
%%Save Depth Data
imgSave=true;

if(imgSave)
     imwrite(OneDresult,[resultFolder ,'Depth\' , name ,'_' , int2str(bands), '-bands_Depth.png']);
     imwrite( ind2rgb(im2uint8(OneDresult), parula(256)), [resultFolder ,'Depth\' , name ,'_' , int2str(bands), '-bands_Depth-imagesc.png']);
end



%%
%%%%% Hierchial cvlustring%%%
addpath('.\lib\hierClustWithLhalf');
cube=mat;%mat(:,:,3:4:end-4);%mat;%datacube;%normat;
cube=normat;
[h,w,matBands]=size(cube);
sz=size(cube);
[h, w, c]=size(cube);

% Select a subset of wavelengths (17) to reduce computational time
% (Remove this line to have the same test as in the paper)
%R = R(:,:,1:10:end); 

% Run H2NMF
algo = 1; 

r = 4; % 4 for others

%[IDX, C, J, sol] = hierclust2nmf(R,r,algo);
[IDX, C, J, sol] = hierclust2nmf(cube,r,algo);


%for i = 1 : size(R,3);
% for i = 1 : size(cube,3);
%     %M(i,:) = reshape(R(:,:,i),sz(1)*sz(2),1);
%     M(i,:) = reshape(cube(:,:,i),sz(1)*sz(2),1);
% end 
% H = nnlsm_blockpivot(C,M); 


M=reshape(cube,[h*w , c]);

H = nnlsm_blockpivot(C,M'); 

nendmembers=reshape(H',[h, w , r]);
%%
[filepath,name,ext] = fileparts(D(k).name);
 for i=1:r
    figure,imagesc(nendmembers(:,:,i));
    %imwrite(normaliseMap(nendmembers(:,:,i)),[name '_Ab_' int2str(i) '.png']);
 end


%%
% saveData=false;
% 
% 
% if(saveData)
%     [filepath,name,ext] = fileparts(curDir(1).name);    
%     save(['results\unmixing\' name '_Data.mat'],'nendmembers','datacube','r');
% end
% for i=1:size(nendmembers,3)
% 
%     figure,imshow(nendmembers(:,:,i),[]);
%     if(saveData)
%         [filepath,name,ext] = fileparts(curDir(i).name);
%         imwrite(normaliseMap(nendmembers(:,:,i)),['results\unmixing\' name '-Abdunce.png']);
%     end
%     
% end
% 
%%

saveImg=true;

if(saveImg)
    for i=1:size(nendmembers,3)

            figure,imshow(nendmembers(:,:,i),[]);

            %[filepath,name,ext] = fileparts(curDir(i).name);
            imwrite(normaliseMap(nendmembers(:,:,i)),[resultFolder ,'Material\' , name ,'_' , int2str(matBands),  int2str(i) '-Abdunce.png']);


    end
end

% [TH,IH]=max(H);
% testimg=reshape(IH,sz(1),sz(2));
% unmixmap=testimg;
% figure,imagesc(unmixmap);

