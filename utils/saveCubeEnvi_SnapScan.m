function [  finalMat, bandName ] = saveCubeEnvi_SnapScan(  filenamepath,dpath,imageflag,savepath )
%SAVECUBEENVI_SNAPSCAN Summary of this function goes here
%   Detailed explanation goes here

if(filenamepath(:,end)=='\')
   
        kind = strfind(filenamepath, '\');
        fname=filenamepath((kind(end-1)+1):end-5);
        opath=filenamepath(1:(kind(end-1)-1));

    else
        kind = strfind(filenamepath, '\');
        fname=filenamepath((kind(end)+1):end-4);
        opath=filenamepath(1:(kind(end)-1));
end

if (nargin > 3)
    
    savefilename=[savepath '\' fname '.mat'];
    if(imageflag)
    imgdir=[savepath '\' fname '_hyperimages\'];
    end
elseif (nargin > 2)
    

    savefilename=[opath '\' fname '.mat'];
    if(imageflag)
    imgdir=[opath '\' fname '_hyperimages\'];
    end

else
    
    
    savefilename=[opath '\' fname '.mat'];
    imageflag=false;
    
    
end

if(imageflag)
    if (exist(imgdir,'dir')==7)
        
    else
        mkdir(imgdir);
    end
end
    


%fname
%savefilename


len = length(filenamepath);
if filenamepath(len-3) =='.'
   prefix = filenamepath(1:len-4);
end
headername = strcat(prefix,'.hdr');
dataname = strcat(prefix,'.raw');

 info = read_envihdr(headername);
 dataC = multibandread(dataname, info.size, [info.format '=>double'], info.header_offset, info.interleave, info.machine);
 disp(['processing ',dataname])
 

norflag=false;

if(info.bands>150)
    norflag=true;%irradiance Needs normalization
    dataCube = double(zeros(info.lines,info.samples,150));
    bandName=strings([150,1]);
    cubeCount = 1;
    for k = 1:info.bands
        if ~isequal(info.wavelength{k},0)
             dataCube(:,:,cubeCount) = dataC(:,:,k);
             bandName(cubeCount)=string(info.band_names{k});
             cubeCount = cubeCount + 1;
        end
    end
else
    dataCube=dataC;
    bandName=info.band_names;
end


%%% Normalization

%%%Dark normalization
if(~isempty(dpath))
    
    %dark_cube=    %not done
    
end
if(norflag)

 dataCube = normalise(dataCube, 'percent', 0.999);
 [numRows,numCols,numSpectra] = size(dataCube);
 scfact = mean(reshape(sqrt(sum(dataCube.^2,3)),numRows*numCols,1));
 dataCube = dataCube./scfact;
 
end


%bandName(2)
if(iscell(bandName))
bandName=cell2mat(bandName);
end    
num = length(bandName);

if(imageflag)
    normat=zeros(size(dataCube));
     for b = 1:num
        slice = dataCube(:,:,b);
        %slice = imadjust(slice);   
        %normat(:,:,b)=slice;


        wl = bandName(b);
        imgname=[fname '_'  num2str(wl) '.png'];
        imgname
        imwrite(slice,[imgdir  '\' imgname]);    
        disp(b);

     end
 
end
 
%%%%%% Reducing Bands
%  finalMat=dataCube(:,:,3:4:end-4);
%  bandName= bandName(3:4:end-4);

figure('Name','Original without rotation or crop','NumberTitle','off'),imshow(dataCube(:,:,12));

finalMat=dataCube;%without roatation
%finalMat=imrotate(dataCube,90); % rotate image

%finalMat=finalMat(849:end-30,:,:); % crop the top part

%figure('Name','After rotation or crop','NumberTitle','off'),imshow(finalMat(:,:,12));

save(savefilename,'finalMat','bandName','-v7.3'); % use this if want to
%save all bands that is 150

%save(savefilename,'finalMat','bandName');

end





