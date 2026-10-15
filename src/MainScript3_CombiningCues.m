% for i=1:r
% 
% figure,imagesc(nendmembers(:,:,i));
% 
% end
[~,fname,fext]=fileparts(D(k).name);
saveFiletemp=[resultFolder fname]
%%

 [ Lap ] = graphLapL2( normat );
 %[ Lap ] = graphLapL2( nendmembers );

%%
finalMats=zeros(size(nendmembers));
i=4;
for i=1:r
% [ Lap2 ] = graphLapL2( nendmembers(:,:,i) );
% fLap=Lap+Lap2;
finalMats(:,:,i) =graphSolveLap( Lap,nendmembers(:,:,i),meanvalpat,nneidis );
figure,imagesc(finalMats(:,:,i));

end

%% not working well

%  allendmembersres=graphSolveLap( Lap,nendmembers,meanvalpat,nneidis );
%  figure,imagesc(allendmembersres);

%%
norFinalMats=zeros(size(nendmembers));
temp=[];
temp(:,:,1)=finalMats(:,:,3);
temp(:,:,2)=finalMats(:,:,4);
temp(:,:,3)=finalMats(:,:,1);

for i=1:r
%if(i~=3)    
%norFinalMats(:,:,i)=normaliseMap(temp(:,:,i)).*i;%+(i*100)
norFinalMats(:,:,i)=normaliseMap(finalMats(:,:,i)).*i;%+(i*100)
%figure,imagesc(norFinalMats(:,:,i));
%end
i
end

res=max(norFinalMats,[],3);
%res=sum(norFinalMats,3);
figure,imagesc(res);


%% Saving files
% for i=1:r
% imwrite( ind2rgb(im2uint8(normaliseMap(finalMats(:,:,i))), parula(256)),  [saveFiletemp '_AbbMaps_' int2str(i) '_imagesc.png']);
% end
testrgbimg=[];
testrgbimg(:,:,1)=normaliseMap(temp(:,:,1));
testrgbimg(:,:,2)=normaliseMap(temp(:,:,2));
testrgbimg(:,:,3)=zeros(size(temp(:,:,2)));%normaliseMap(temp(:,:,3));%zeros(size(temp(:,:,2)));%zeros(size(nendmembers(:,:,1)));%norFinalMats(:,:,3);
figure,imshow(testrgbimg);
% 
% imwrite( ind2rgb(im2uint8(normaliseMap(res)), parula(256)),  [saveFiletemp '_final_imagesc.png']);
 imwrite( testrgbimg,  [saveFiletemp '_final_RGB.png']);

%% drawing bounding box

normat=datacube;%need to be commented
test=[];

test(:,:,1)=normaliseMap(temp(:,:,1));
test(:,:,2)=normaliseMap(temp(:,:,2));
%test(:,:,3)=normaliseMap(temp(:,:,3));
% [val,mind]=max(test(:));
% val
% [ht, wt]=size(test);
% [row,col] = ind2sub([ht, wt],mind);
% T = graythresh(I);
% T
% for i=10:10:500
%     i
%     cmaxval=test(col+i,row)
%     if(cmaxval<(T-0.05))
%         cmax=col+i
%         break;
%     end
% end
eColor=['r','g','b','c','m','y','k','w'];
figure,imshow(normat(:,:,30));

for i=1:size(test,3)%r-1
    BW1 = imbinarize(test(:,:,i));
    BW2 = bwareafilt(BW1, 1);
    stats = regionprops(BW2,'BoundingBox');
    recdim=stats.BoundingBox;
    %rectangle('Position',[recdim(1),recdim(2),recdim(3)+20,recdim(4)+20],'EdgeColor',eColor(i),'LineWidth',2);
    rectangle('Position',recdim,'EdgeColor',eColor(i),'LineWidth',2);
end

%figure,imshow(finalMats(:,:,1));
% figure,imshow(normat(:,:,12));
% rectangle('Position', [673-200,840-120,350,350],'EdgeColor','r','LineWidth',2);
