%% STIPS algorithm enables tracking labyrinthine patterns and reveals distinct rhythmic dynamics of actin microridges
%% Authors: Rajasekaran Bhavna1*, Mahendra Sonawane1
%% 1 Department of Biological Sciences, Tata Institute of Fundamental Research, Colaba, Mumbai- 400005

function [SAL,t1_in, t2_in,stats_t1, stats_t2,dbscanpnts,dbradius,creatradius]= create_intermediate_temporal_links(stats_t1, stats_t2,t2_joins_binar)
clearvars -except stats_t1 stats_t2 t1_binar t2_joins_binar 
dbscanpnts=1;
dbradius=1;
creatradius=2;

%% find postions of zeros at t2
[pos0(:,1),pos0(:,2)]=find(t2_joins_binar==0);

%% find postions of ones at t2
[pos1(:,1),pos1(:,2)]=find(t2_joins_binar==1);

clear stats_break_t1;
for e=1:length(stats_t1)
    stats_break_t1(e).PixelList=stats_t1(e).PixelList;
end

%% set pos0 t2 postions to zero at t1
for s=1:length(stats_break_t1)
    dd=find(ismember(stats_break_t1(s).PixelList,pos0,'rows'));
        for y=1:length(dd)
            stats_break_t1(s).PixelList(dd(y),:)=[NaN NaN];
        end
        stats_break_t1(s).PixelList=rmmissing(stats_break_t1(s).PixelList,1);  
    clear dd;
end
clear CPL;
CPL=struct2cell(stats_break_t1);
CPL=squeeze(CPL);
CPL=transpose(CPL);
distparam='euclidean';


%% to add pos1 t2 positions to t1, using closest distance
clear valind minvalind
for m=1:length(pos1)
    clear X MdlES 
    X=pos1(m,:);
    MdlES = createns(X,'NSMethod','exhaustive');
    %MdlES = createns(X,'NSMethod','kdtree','BucketSize',200);
   clear IdxNN D IdNN DxNN med_dist
    for n=1:length(CPL)
        [IdxNN{n},D{n}] = rangesearch(MdlES,CPL{n},creatradius);
        IdNN{n}=IdxNN{n}(~cellfun('isempty',IdxNN{n}));
        DxNN{n}=D{n}(~cellfun('isempty',D{n}));
        med_dist{n}=mean(cell2mat(DxNN{n}),'omitnan');
    end
    valind{m}=find(~cellfun(@isempty,IdNN));
    [minvalind(m,1),minvalind(m,2)]=min(cell2mat(med_dist));
end
%% add all unassigned_pos1 at the end of cell array
clear unassigned_pos1 unassignedpos idx_unassign max_unassign

for m=1:length(pos1)
    if isnan(minvalind(m,1))
        unassigned_pos1(m,:)=pos1(m,:);
    else
        CPL{minvalind(m,2)}=cat(1,CPL{minvalind(m,2)},pos1(m,:));
    end
end
unassignedpos = unassigned_pos1(any(unassigned_pos1,2),:);

%% check within unassignedpos for interdist matches

idx_unassign = dbscan(unassignedpos,dbradius,dbscanpnts,'Distance',distparam);
max_unassign=max(idx_unassign);
clust_unassign={};
for ui=1:max_unassign
    clust_unassign{end+1}=unassignedpos(idx_unassign==ui,:);
end

clear PPL
for g=1:length(CPL)
    PPL{g}=sortrows(CPL{g});
end

%% split a pixcluster within the List
IPL={};
for g=1:length(PPL)
    clear Y idx mimax xx
    if ~isempty(PPL{g})
        Y=PPL{g};
        idx = dbscan(Y,dbradius,dbscanpnts,'Distance',distparam);
        mimax=max(idx);
        if mimax>1
            for mi=2:mimax
                IPL{end+1}=PPL{g}(idx==mi,:);
            end
        end
        xx=find(idx~=1);
        for x=1:length(xx)
            PPL{g}(xx(x),:)=NaN;
        end
    end
end

for k=1:length(PPL)
PPL{k}=rmmissing(PPL{k});
PPL{k}= unique(PPL{k},'rows','stable');
end

for h=1:length(IPL)
    IPL{h}= unique(IPL{h},'rows','stable');
end

for k=1:length(clust_unassign)
    clust_unassign{k}=unique(clust_unassign{k},'rows','stable');
end
clear newCL SAL
newCL=cat(2,PPL,IPL,clust_unassign);
newCL=newCL(~cellfun('isempty',newCL));
SAL = cell2struct(newCL,'PixelList',length(newCL));


for r=1:length(SAL)
    SAL(r).PixelList=flip(SAL(r).PixelList,2);
end


%% consolidate SAL as per stats_t1
for ls1=1:length(SAL)
    for lt1=1:length(stats_t1)
        t1_in{ls1,lt1}=find(ismember(stats_t1(lt1).PixelList,SAL(ls1).PixelList,'rows'));
    end
end

%% consolidate SAL as per stats_t2
for ls2=1:length(SAL)
    for lt2=1:length(stats_t2)
        t2_in{ls2,lt2}=find(ismember(stats_t2(lt2).PixelList,SAL(ls2).PixelList,'rows'));
    end
end

[r1,c1]=size(t1_in);
[r2,c2]=size(t2_in);

if r2 ~= r1
    n=abs(r2-r1);
    if r1>r2
        t2_in=cat(1,t2_in,cell(n,c2));
    else
        t1_in=cat(1,t1_in,cell(n,c1));
    end
end

%% plots
%{
figure;
for l1=1:length(stats_t1)
    plot(stats_t1(l1).PixelList(:,1),stats_t1(l1).PixelList(:,2),'g*','LineWidth',2);
    hold on;
    text(stats_t1(l1).PixelList(1,1),stats_t1(l1).PixelList(1,2),num2str(l1),'FontSize', 18,'Color','blue');
    hold on;
end
title('t1')
gca.XDir = 'reverse';
gca.YDir = 'reverse';
view(0,-90);
box off;
axis off;
hold off;

figure;
for l1=1:length(SAL)
    plot(SAL(l1).PixelList(:,1),SAL(l1).PixelList(:,2),'m*','LineWidth',2);
    hold on;
    x=find(~isnan(SAL(l1).PixelList(:,1)));
    if ~isempty(x)
        text(SAL(l1).PixelList(x(1),1),SAL(l1).PixelList(x(1),2),num2str(l1),'FontSize', 15);
    end
    hold on;
end
title('intermediate t1-t2')
gca.XDir = 'reverse';
gca.YDir = 'reverse';
view(0,-90);
 daspect([1 1 1])
box off;
axis off;
hold off;
%clear gca;
%clear f1;
%f1=gca;exportgraphics(f1,['SAL',num2str(radiuspix),'_',num2str(dbscanpnts),'.eps'],'Resolution',600);close;clear gca;

figure;
for l2=1:length(stats_t2)
    plot(stats_t2(l2).PixelList(:,1),stats_t2(l2).PixelList(:,2),'c*','LineWidth',2);
    hold on;
    text(stats_t2(l2).PixelList(1,1),stats_t2(l2).PixelList(1,2),num2str(l2),'FontSize', 18,'Color','black');
    hold on;
end
title('t2')
gca.XDir = 'reverse';
gca.YDir = 'reverse';
view(0,-90);
box off;
axis off;
hold off;
%}
