%% STIPS algorithm enables tracking labyrinthine patterns and reveals distinct rhythmic dynamics of actin microridges
%% Authors: Rajasekaran Bhavna1*, Mahendra Sonawane1
%% 1 Department of Biological Sciences, Tata Institute of Fundamental Research, Colaba, Mumbai- 400005

function cellimg = readpng_rgb(imgpng)
warning('off');
filenamepath=regexp(imgpng,'/','split');
filename=regexp(filenamepath{end},'\.','split');
if strcmp(filename{2},'png')
    cellimg=imread(imgpng);
    Im = imresize(cellimg,[256 256]);
    Irep = repmat(Im,[1 1 1]);
    cellimg=Irep;
end
