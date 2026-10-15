

%DCube=diffcube;%GRcube;
DCube=mat;%gedgecube;%divcube;%GRcube;%testmat;%TM;%diffDcube;%diffSapGC;
%DCube=diffGMZ;


f = figure;
hn = uicontrol('Position',[20 20 200 40],'String','next',...
              'Callback','uiresume(gcbf)');
hp = uicontrol('Position',[230 20 200 40],'String','previous',...
              'Callback','uiresume(gcbf)');

 htxt = uicontrol('Style','text',...
        'Position',[460 20 200 25],...
        'String','Dummy Value');
    
          
  sz=size(DCube);        
%disp('This will print immediately');


 i=1;
 %htxt.String=bandName(i);
 set(htxt,'String',['',bandName(i)]);
while true
    
  

   %imagesc(DCube(:,:,i));
   imshow(DCube(:,:,i));
   %pause(6);
   
   uiwait(gcf); 
   
   
   
   %disp('This will print after you click Continue');
  
   if ~ishandle(f)
       i=sz(3)+1;
       break;
   end
   
   han = gco(f);
    if han==hn
       i=i+1;
    else
        i=i-1;
    end
    
   if i==sz(3)+1
   i=1;
   elseif i==0
   i=sz(3);
   end
   
    set(htxt,'String',['',bandName(i)]);
   i
end