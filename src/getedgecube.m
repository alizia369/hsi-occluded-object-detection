function [edgesum,edgecube]=getedgecube(normat)

%GETEDGECUBE Summary of this function goes here
%   Detailed explanation goes here

sz=size(normat);

edgecube=zeros(sz);
normat=removenoise(normat);

for i = 1:sz(3)
    
    slice = normat(:,:,i);
    %slice = imadjust(slice);   
    
    edgeMap=edge(slice,'canny',0.30);%previous 0.25

    edgecube(:,:,i)=edgeMap;
     
%     disp('edge');
%     i
end

edgesum=sum(edgecube,3);

%edgesum(edgesum>1)=1;

%figure,imshow(sumedg);




end

