addpath('.\lib\hierClustWithLhalf');
cube=mat(:,:,3:4:end-4);%datacube;%normat;
sz=size(cube);
[h, w, c]=size(cube);

% Select a subset of wavelengths (17) to reduce computational time
% (Remove this line to have the same test as in the paper)
%R = R(:,:,1:10:end); 

% Run H2NMF
algo = 1; r = 6; % 4 for others
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
saveData=false;


if(saveData)
    [filepath,name,ext] = fileparts(curDir(1).name);    
    save(['results\unmixing\' name '_Data.mat'],'nendmembers','datacube','r');
end
for i=1:size(nendmembers,3)

    figure,imshow(nendmembers(:,:,i),[]);
    if(saveData)
        [filepath,name,ext] = fileparts(curDir(i).name);
        imwrite(normaliseMap(nendmembers(:,:,i)),['results\unmixing\' name '-Abdunce.png']);
    end
    
end

%%

saveImg=true;

if(saveImg)
    for i=1:size(nendmembers,3)

        %figure,imshow(nendmembers(:,:,i),[]);

            %[filepath,name,ext] = fileparts(curDir(i).name);
            imwrite(normaliseMap(nendmembers(:,:,i)),['.\data\snapscan b150u\Results\Material\cube3_36bands_' int2str(i) '-Abdunce.png']);


    end
end

% [TH,IH]=max(H);
% testimg=reshape(IH,sz(1),sz(2));
% unmixmap=testimg;
% figure,imagesc(unmixmap);

%% Testing

% endmembers(:,:,1)=reshape(H(1,:),sz(1),sz(2));
% endmembers(:,:,2)=reshape(H(2,:),sz(1),sz(2));
% endmembers(:,:,3)=reshape(H(3,:),sz(1),sz(2));
% endmembers(:,:,4)=reshape(H(4,:),sz(1),sz(2));
% endmembers(:,:,5)=reshape(H(5,:),sz(1),sz(2));
% endmembers(:,:,6)=reshape(H(6,:),sz(1),sz(2));
% endmembers(:,:,7)=reshape(H(7,:),sz(1),sz(2));
% endmembers(:,:,8)=reshape(H(8,:),sz(1),sz(2));
% endmembers(:,:,9)=reshape(H(9,:),sz(1),sz(2));
% 
% for i=1:9
% figure,imshow(endmembers(:,:,i),[]);
% 
% end