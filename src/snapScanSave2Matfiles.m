mainFolder='.\data\snapscan b150u\6\';
dirFolder=[mainFolder,'*.hdr'];
D = dir(dirFolder); % A is a struct ... first elements are '.' and '..' used for navigation.
 
for i = 1:length(D) % print folders and k values
 %k   
 %D(k).name  %Get the current subdirectory name
  [  finalMat, bandName ]=saveCubeEnvi_SnapScan([mainFolder , D(i).name],'' );
   %save([mainFolder,D(i).name,'.mat'],'finalMat','bandName','-v7.3');%already_saved_in_saveCubeEnvi_SnapScan(---)
   imwrite(finalMat(:,:,30),[mainFolder , D(i).name , '_30band.png']);

  disp([num2str(i) '. ' D(i).name]);
end

%% test


% [  finalMat, bandName ]=saveCubeEnvi_SnapScan([mainFolder , D(i).name],'' );
%    %save([mainFolder,D(i).name,'.mat'],'finalMat','bandName','-v7.3');%already_saved_in_saveCubeEnvi_SnapScan(---)
%    imwrite(finalMat(:,:,30),[mainFolder , D(i).name , '_30band.png']);
%    disp([num2str(i) '. ' D(i).name]);


%%

% finalMat=finalMatChk1;
% [h,w,c]=size(finalMat);
% for i = 1:c % print folders and k values
%  
%     finalMat(:,:,i)=imadjust(finalMatChk1(:,:,i));
% end
% 
% figure,imshow(finalMat(:,:,12));
% 
% save('.\data\snapscan b150u\3ourside\outdoor1_corrected.mat','finalMat','bandName','-v7.3');

%%

% mainFolder='.\data\snapscan b150u\3ourside\';
% resultFolder='.\data\snapscan b150u\3ourside\Results\'
% D = dir([mainFolder , '*.mat']); % A is a struct ... first elements are '.' and '..' used for navigation.
%  
% for i = 1:length(D) % print folders and k values
%  
%  disp([num2str(i) '. ' D(i).name]);
%  name=D(i).name;
%  load([mainFolder , name]);
%  imwrite(finalMat(:,:,30),[mainFolder , name , '_30band.png']);
%  
% end

