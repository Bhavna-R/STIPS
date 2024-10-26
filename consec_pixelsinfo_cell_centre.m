%% STIPS algorithm enables tracking labyrinthine patterns and reveals distinct rhythmic dynamics of actin microridges
%% Authors: Rajasekaran Bhavna1*, Mahendra Sonawane1
%% 1 Department of Biological Sciences, Tata Institute of Fundamental Research, Colaba, Mumbai- 400005

function [stats_t1, stats_t2,t2_joins_binar,t2_joinsbyt1]= consec_pixelsinfo_cell_centre(timeframe,Rf1,px_xysz)
clearvars -except timeframe Rf1 px_xysz

timeframenext=timeframe+1;

cell_centre=1;
% for whole cell, set cell_centre=0;
if cell_centre==1 % cell cropped region
  t1=imcrop((Rf1{1, timeframe}),[75 50 100 100]); 
  t2=imcrop((Rf1{1, timeframenext}),[75 50 100 100]);
else % whole cell
  t1=Rf1{1, timeframe};
  t2=Rf1{1, timeframenext};
end

t2byt1=t2./t1;

t2_replace=t2;

t2byt1_replace=t2byt1;

t2_replace(t2_replace==0)=NaN;
logt2_replace=log(t2_replace);

t2_joins=t2byt1_replace+(1-logt2_replace);
t2_joins(isnan(t2_joins))=0;
t2_joins(isinf(t2_joins))=5;

t1_binar=imbinarize(t1);
t2_binar=imbinarize(t2);
t2_joins_binar=imbinarize(t2_joins);

%0 breaks,inf joins, NaNs nochange in 0's, 1's nochange in common pix
t2_joinsbyt1=t2_joins_binar./t1_binar;
[t2_jns(:,1),t2_jns(:,2)]=find(t2_joinsbyt1==inf);
[t2_brks(:,1),t2_brks(:,2)]=find(t2_joinsbyt1==0);
[t1_t2_nochng0(:,1),t1_t2_nochng0(:,2)]=find(isnan(t2_joinsbyt1)==1);
[t1_t2_comn1pix(:,1),t1_t2_comn1pix(:,2)]=find(t2_joinsbyt1==1);

BW_t1 = bwareaopen(t1_binar,1);
BW_t2=bwareaopen(t2_joins_binar,1);

CC1 = bwconncomp(BW_t1);
CC2 = bwconncomp(BW_t2);

stats_t1=regionprops(CC1,'PixelList');%pix at t1
stats_t2=regionprops(CC2,'PixelList');%pix at t2
