%% STIPS algorithm enables tracking labyrinthine patterns and reveals distinct rhythmic dynamics of actin microridges
%% Authors: Rajasekaran Bhavna1*, Mahendra Sonawane1
%% 1 Department of Biological Sciences, Tata Institute of Fundamental Research, Colaba, Mumbai- 400005

function [inprsz,pxdsResults,px_xysz]=rescalepixsiz_temporal(locateraw,nw)
pixmic=0.1977;
dirraw = fullfile(locateraw,'*.png');
fileraw = dir(dirraw);

clear rwimg padarr
 for f=1:length(fileraw)
     rwimg{f}=imread([fileraw(f).folder,'/',fileraw(f).name]);
end
[rsz csz]=cellfun(@size,rwimg,'UniformOutput',false);
rz=max(cell2mat((rsz)));
cz=max(cell2mat((csz)));
maxrc=max(rz,cz);
if (mod(maxrc,2)~=0),maxrc=maxrc+1;end
topadarray=rwimg;

for f=1:length(topadarray)
    
    if (mod(rsz{f},2)~=0 && mod(csz{f},2)==0)
        zz=zeros(1,csz{f});
        topadarray{f}=cat(1,topadarray{f},zz);
        
    elseif (mod(rsz{f},2)==0 && mod(csz{f},2)~=0)
        zz=zeros(rsz{f},1);
        topadarray{f}=cat(2,topadarray{f},zz);
        
    elseif (mod(rsz{f},2)~=0 && mod(csz{f},2)~=0)
        zz=zeros(1,csz{f});
        topadarray{f}=cat(1,topadarray{f},zz);
        [nr nc]=size(topadarray{f});
        zz=zeros(nr,1);
        topadarray{f}=cat(2,topadarray{f},zz);
    end
    
    [pr pc]=size(topadarray{f});
    add_r=(maxrc-pr)/2;
    add_c=(maxrc-pc)/2;
    padarr{f}=padarray(topadarray{f},[add_r add_c]); 
    padsz=size(padarr{f});
    inprsz{f}=imresize(padarr{f},[256 256]);
    pxdsResults{f} = semanticseg(inprsz{f},nw,'OutputType','categorical');
end
px_xysz(1)=(padsz(1)*pixmic)./256;
