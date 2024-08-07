%% STIPS algorithm enables tracking labyrinthine patterns and reveals distinct rhythmic dynamics of actin microridges
%% Authors: Rajasekaran Bhavna1*, Mahendra Sonawane1
%% 1 Department of Biological Sciences, Tata Institute of Fundamental Research, Colaba, Mumbai- 400005

clear all;
userpath('reset');

load('comb14_nw93.mat');
nw=Nw14_trainedNetworks93;
layer='Softmax-Layer';
locateraw='Raw_Icell3_011019yolkcentre52hpf';

    time_inter_sec=0.56; %(seconds)
    datasetdir=locateraw;
    imagedir=fullfile(datasetdir);
    imdsseries = imageDatastore(imagedir);
    Imageseries = imdsseries.Files;
    imdsseries = imageDatastore(Imageseries,'FileExtensions','.png','ReadFcn',@readpng_rgb);
    imgseries = imdsseries.Files;
    [~,pxdsResults,px_xysz]=rescalepixsiz_temporal(locateraw,nw);

    featuresseries = activations(nw,imdsseries,layer);
    szf=size(featuresseries);%256,256,f=2,t

    for t=1:szf(4)
        FT1(:,:,t) = featuresseries(:,:,1,t);
        FT2(:,:,t) = featuresseries(:,:,2,t);
    end
    ft1_th=FT1.*(FT1>0.5);
    ft2_th=FT2.*(FT2>0.5);

    for t=1:szf(4)
        Rf1{t}=ft1_th(:,:,t);
        Rf2{t}=ft2_th(:,:,t);
    end

distparam='euclidean';
radiuspix=3;

for timeframe=1:(szf(4)-1)

    disp(timeframe);

    clearvars -except timeframe timeframenext Rf1 time_inter_sec pxdsResults px_xysz Link_list time_inter_sec micron_scale szf creatradius radiuspix a distparam
    clear stats_t1 stats_t2 t2_joins_binar SAL t1_in t2_in new_stats_asgn Assign_tbl A_srttbl two_frm_lnks
    [stats_t1,stats_t2,t2_joins_binar,t2_joinsbyt1]= consec_pixelsinfo_cell_centre(timeframe,Rf1,px_xysz);
    
    [SAL,t1_in,t2_in,stats_t1,stats_t2,dbscanpnts,dbradius,creatradius]= create_intermediate_spatiotemporal_links(stats_t1,stats_t2,t2_joins_binar);

    clear f1 f2 fr1 fc1 fr2 fc2 B C disappear_t1
    f1=(~cellfun('isempty',t1_in));
    f2=(~cellfun('isempty',t2_in));
    [fr1,fc1]=size(f1);
    [fr2,fc2]=size(f2);

    B=cellfun(@any,t1_in);
    C=any(B,1);
    disappear_t1=find(C==0);

    clear tprev  tnex tpn
    for fr=1:fr1
        tprev{fr,1}=find(f1(fr,:));
        tnex{fr,1}=find(f2(fr,:));

        tpn{fr,1}=find(f1(fr,:));
        tpn{fr,2}=find(f2(fr,:));
    end

    % check if entire row is 0 & remove such rows
    As=cellfun(@isempty,tpn);
    TF=find(all(As == 1,2));
    tpn(TF,:) = [];
    tprev(TF,:) = [];
    tnex(TF,:) = [];

    new_stats_asgn=[];
    saml_don_id=[];
    ind_m=[];
    appear_ind=find(cellfun(@isempty,tprev));
    y_ind=1;

    for y=1:length(tpn)
        if (ismember(y,saml_don_id)==1)
            continue;
        end
        % disp(y);
        % disp('**');
        % disp(length(saml_don_id));
        % disp('**');

        tp=tpn{y,1};
        tn=tpn{y,2};
        clear m2tp m2tn ind_m len_m2tn;
        unq_m2tp_tn=[];
        unq_m2tn_tp=[];
        appear_ind=setxor(appear_ind,intersect(saml_don_id,appear_ind));

        m2tp={};
        m2tn={};
        for r=1:length(tp)
            m2tp{r}=find(cell2mat(cellfun(@(x) any(x==tp(r)), tprev,'UniformOutput',false)));
        end
        for r=1:length(tn)
            m2tn{r}=find(cell2mat(cellfun(@(x) any(x==tn(r)), tnex,'UniformOutput',false)));
        end
        clear m2tp_m2tn
        m2tp_m2tn=unique(cat(1,m2tp{:},m2tn{:}));

        %% for the existing m2tp & m2tn list check for other potential matches with tpn array
        if (exist('m2tp_m2tn') & ~isempty(m2tp_m2tn))
            clear pot_matches_prev pot_matches_nex
            for r=1:length(m2tp_m2tn)
                for zz=1:length(cell2mat(tpn(m2tp_m2tn(r),1)))
                    clear z;
                    z=cell2mat(tpn(m2tp_m2tn(r),1));
                    pot_matches_prev{r}{zz}=find(cell2mat(cellfun(@(x) any(x==z(zz)), tprev,'UniformOutput',false)));
                end
                for zz=1:length(cell2mat(tpn(m2tp_m2tn(r),2)))
                    clear z;
                    z=cell2mat(tpn(m2tp_m2tn(r),2));
                    pot_matches_nex{r}{zz}=find(cell2mat(cellfun(@(x) any(x==z(zz)), tnex,'UniformOutput',false)));
                end
            end

            clear m2tp_prev
            if (exist('pot_matches_prev'))
                for s=1:length(pot_matches_prev)
                    if ~isempty(pot_matches_prev{s})
                        m2tp_prev{s}=cat(1,pot_matches_prev{s}{:});
                    end
                end
                m2tp_prev=unique(cat(1,m2tp_prev{:}));
                clear m2tp
                m2tp=m2tp_prev;
            end

            clear m2tn_nex
            if (exist('pot_matches_nex'))
                for s=1:length(pot_matches_nex)
                    if ~isempty(pot_matches_nex{s})
                        m2tn_nex{s}=cat(1,pot_matches_nex{s}{:});
                    end
                end

                m2tn_nex=unique(cat(1,m2tn_nex{:}));
                clear m2tn
                m2tn=m2tn_nex;
            end

            append_out=[];
            append_out= cat(1,m2tp,m2tn);
            if  iscell(append_out), append_out=cat(1,append_out{:});end
            m2tp_prev_m2tn_nex=unique(append_out);
        end

        m2tp_final=[];
        m2tn_final=[];

        for p=1:20
            clear pot_matches_prev pot_matches_nex
            for r=1:length(m2tp_prev_m2tn_nex)
                for zz=1:length(cell2mat(tpn(m2tp_prev_m2tn_nex(r),1)))
                    clear z;
                    z=cell2mat(tpn(m2tp_prev_m2tn_nex(r),1));
                    pot_matches_prev{r}{zz}=find(cell2mat(cellfun(@(x) any(x==z(zz)), tprev,'UniformOutput',false)));
                end
                for zz=1:length(cell2mat(tpn(m2tp_prev_m2tn_nex(r),2)))
                    clear z;
                    z=cell2mat(tpn(m2tp_prev_m2tn_nex(r),2));
                    pot_matches_nex{r}{zz}=find(cell2mat(cellfun(@(x) any(x==z(zz)), tnex,'UniformOutput',false)));
                end
            end
            clear m2tp_prev m2tn_nex
            if exist('pot_matches_prev')
                for s=1:length(pot_matches_prev)
                    if ~isempty(pot_matches_prev{s})
                        m2tp_prev{s}=cat(1,pot_matches_prev{s}{:});
                    else
                        m2tp_prev{s}=  m2tp_prev_m2tn_nex(s);
                    end
                end
            else
                m2tp_prev=m2tp_prev_m2tn_nex;
            end

            if exist('pot_matches_nex')
                for s=1:length(pot_matches_nex)
                    if ~isempty(pot_matches_nex{s})
                        m2tn_nex{s}=cat(1,pot_matches_nex{s}{:});
                    else
                        m2tn_nex{s}=m2tp_prev_m2tn_nex(s);
                    end
                end
            else
                m2tn_nex=m2tp_prev_m2tn_nex;
            end

            if  iscell(m2tp_prev)
                m2tp_prev=unique(cat(1,m2tp_prev{:}));
            end
            m2tp_prev_m2tn_nex=unique(cat(1,m2tp_prev_m2tn_nex,m2tp_prev));
            if exist('m2tp_final') & ~isempty(m2tp_final)
                m2tp_final=unique(cat(1,m2tp_final,m2tp_prev_m2tn_nex));
            else
                m2tp_final =m2tp_prev_m2tn_nex;
            end
            if  iscell(m2tn_nex)
                m2tn_nex=unique(cat(1,m2tn_nex{:}));
            end
            m2tp_prev_m2tn_nex=unique(cat(1,m2tp_prev_m2tn_nex,m2tn_nex));
            if exist('m2tn_final') &~isempty(m2tn_final)
                m2tn_final=unique(cat(1,m2tn_final,m2tp_prev_m2tn_nex));
            else
                m2tn_final =m2tp_prev_m2tn_nex;
            end
            continue;
        end
        m2tp_final=unique(m2tp_final);
        m2tn_final=unique(m2tn_final);

        m2tp_m2tn_final=[];
        if  iscell(m2tp_final), m2tp_final=cat(1,m2tp_final{:});end
        if  iscell(m2tn_final), m2tn_final=cat(1,m2tn_final{:});end
        m2tp_m2tn_final=unique(cat(1,m2tp_final,m2tn_final));

        clear new_matches unqnew_matches

        %% to addon matches from m2tp to tnex
        if (exist('m2tp_final') & ~isempty(m2tn_final))
            for e=1:length(m2tp_m2tn_final)
                clear z
                z=cell2mat(tnex(m2tp_m2tn_final(e,1)));
                for zz=1:length(z)
                    new_matches{zz}{e}=find(cell2mat(cellfun(@(x) any(x==z(zz)), tnex,'UniformOutput',false)));
                end
            end
            if (exist('new_matches') & ~isempty(new_matches))
                for g=1:length(new_matches),new_matches{g}=cat(1,new_matches{g}{:});end
                unqnew_matches= unique(cat(1,new_matches{:}));
                m2tn_final=unique(cat(1,[m2tn_final],unqnew_matches));
            end
        end

        %%
        clear new_matches unqnew_matches

        %% to addon matches from m2tn to tprev
        if (exist('m2tn_final') & ~isempty(m2tp_final))
            for e=1:length(m2tp_m2tn_final)
                clear z
                z=cell2mat(tprev(m2tp_m2tn_final(e,1)));
                for zz=1:length(z)
                    new_matches{zz}{e}=find(cell2mat(cellfun(@(x) any(x==z(zz)), tprev,'UniformOutput',false)));
                end
            end
            if (exist('new_matches') & ~isempty(new_matches))
                for g=1:length(new_matches),new_matches{g}=cat(1,new_matches{g}{:});end
                unqnew_matches= unique(cat(1,new_matches{:}));
                m2tp_final=unique(cat(1,[m2tp_final],unqnew_matches));
            end
        end

        Finalm2tp_m2tn_final=[];
        Finalm2tp_m2tn_final=unique(cat(1,m2tp_final,m2tn_final));

        check_app=(tnex(appear_ind));
        val_unq_tptn=tnex(Finalm2tp_m2tn_final);

        QR =  1:length(check_app);
        CR= 1:length(val_unq_tptn);

        clear checkapp ind_take consider_appear
        for  app_quer= QR
            for curnt_list=CR
                [checkapp{app_quer,curnt_list}] = find(intersect([check_app{app_quer}],[val_unq_tptn{curnt_list}],'stable'));
            end
        end
        if (exist('checkapp'))
            checkapp= ~cellfun('isempty', checkapp);
            ind_take=find(sum(checkapp,2));%column vector containing the sum of each row
            consider_appear=appear_ind(ind_take);
            concat_appear=[];
            concat_appear=cat(1,concat_appear,consider_appear);
            Finalm2tp_m2tn_final=cat(1,Finalm2tp_m2tn_final,concat_appear);
        end
        clear unq_tptn
        unq_tptn=unique(Finalm2tp_m2tn_final);

        all_tpj=[];
        all_tnb=[];
        clear tpj tnb pix_tpj pix_tnb t1_rc t2_rc;
        clear join_pix_tpj fin_pix_tpj centpx_tp clustsize clustid;
        clear unq_all_tpj  unq_all_tnb len_unq_all_tnb nonemp_t1_rc act_t2_rc;

        for l=1:length(unq_tptn)
            tpj=tpn{unq_tptn(l),1};
            tnb=tpn{unq_tptn(l),2};
            if ~(isempty(tpj))
                for p=1:length(tpj)
                    pix_tpj{l}{p}=stats_t1(tpj(p)).PixelList(t1_in{unq_tptn(l),tpj(p)},:);
                end
            else
                pix_tpj{l}={[]};
            end
            for q=1:length(tnb)
                pix_tnb{l}{q}=stats_t2(tnb(q)).PixelList(t2_in{unq_tptn(l),tnb(q)},:);
            end
        end

        for l=1:length(unq_tptn)
            clear tpj tnb;
            tpj=tpn{unq_tptn(l),1};
            tnb=tpn{unq_tptn(l),2};
            join_pix_tpj{l} = cat(1,pix_tpj{l}{:});

            t1_rc{l}=tpj;
            t2_rc{l}=tnb;

            if (length(t2_rc{l})==1 && ~isempty(join_pix_tpj{l}))
                fin_pix_tpj{l} = cat(1,join_pix_tpj{l});

                if (numel(fin_pix_tpj{l})==2)
                    centpx_tp{l}=fin_pix_tpj{l};
                else
                    centpx_tp{l}=mean(fin_pix_tpj{l});
                end

            elseif (length(t2_rc{l})>1 && ~isempty(join_pix_tpj{l}))
                clear clustsize clustid max_cid dist_temp S;
                clustsize=length(pix_tnb{l});
            

                if (length(join_pix_tpj{l}(:,1)) <clustsize)
                    clear appd rr;
                    appd=clustsize-length(join_pix_tpj{l}(:,1))+1;
                    if mod(appd,2)~= 0
                        appd=appd+1;
                    end
                    rr=repmat(join_pix_tpj{l},appd,1);
                    join_pix_tpj{l}=rr;

                    %clustsize=length(join_pix_tpj{l}(:,1));
                end
                clear xccord ycoord std_x std_y newcoord newradii
                xccord=join_pix_tpj{l}(:,1);
                ycoord=join_pix_tpj{l}(:,2);
                std_x=std(xccord);
                std_y=std(ycoord);

                if ((std_x<=3.5 && std_x>0) && (std_y<=3.5 && std_y>0))
                    newcoord=[xccord ycoord];
                elseif (std_x==0 || std_y==0)
                    newcoord=[xccord ycoord];
                else
                    newcoord=[xccord./std_x ycoord./std_y];
                end
                dist_temp =squareform(pdist(join_pix_tpj{l}));
                newradii=floor(sqrt(mean(mean(dist_temp))));
                if newradii==0
                    newradii=newradii+1;% rare case with only 1 element
                end
                if (newradii<radiuspix)
                    try
                        rng('default');
                         clustid =spectralcluster(newcoord,clustsize,'SimilarityGraph','epsilon','Radius',newradii);

                    catch
                        newradii=newradii+1;
                        rng('default');
                              end
                    clustid =spectralcluster(newcoord,clustsize,'SimilarityGraph','epsilon','Radius',newradii);
                else
                    rng('default');
                            clustid =spectralcluster(newcoord,clustsize,'SimilarityGraph','epsilon','Radius',radiuspix);
                end

                max_cid=max(clustid);
                for mi=1:max_cid
                    fin_pix_tpj{l}{mi}(:,:)=join_pix_tpj{l}(clustid==mi,:);
                    if (numel(fin_pix_tpj{l}{mi})==2)
                        centpx_tp{l}(mi,1:2)=fin_pix_tpj{l}{mi};
                    else
                        centpx_tp{l}(mi,1:2)=mean(fin_pix_tpj{l}{mi});
                    end
                end
           
                clear clustsize clustid max_cid;
                %*******************************************************************
            elseif (length(t2_rc{l})>=1 && isempty(join_pix_tpj{l}))
                fin_pix_tpj{l} =cat(1,join_pix_tpj{l});
                centpx_tp{l}=[];

                %elseif (length(t2_rc{l})==1 && isempty(join_pix_tpj{l}))
                %   fin_pix_tpj{l} =cat(1,join_pix_tpj{l});
                %  centpx_tp{l}=[];
            end
            all_tpj=cat(2,all_tpj,tpj);
            all_tnb=cat(2,all_tnb,tnb);
        end

        for s=1:length(fin_pix_tpj)
            if (iscell(fin_pix_tpj{s})==0 && isempty(fin_pix_tpj{s})),fin_pix_tpj{s}=cell(0,0);end
        end

        unq_all_tpj=unique(all_tpj);
        unq_all_tnb=unique(all_tnb);
        len_unq_all_tnb=length(unq_all_tnb);

        %nonemp_t1_rc=find(~cellfun('isempty',t1_rc));
        %act_t1_rc=t1_rc(nonemp_t1_rc);

        act_t2_rc=t2_rc;
        all_pix_tpj=[];
        joinpix_tnb=[];
        all_pix_tnb=[];

        %% make as many t1_rc as t2_rc using fin_pix_tpj (after spectralcluster re-assign index)
        clear new_t1_rc act_t1_rc
        for kk=1:length(fin_pix_tpj)
            clear ss bb lp s fp b sss bbb
            if (~isempty(fin_pix_tpj{kk}) &&  iscell(fin_pix_tpj{kk}))
                lp=cellfun(@numel,pix_tpj{kk});
                s=find((lp<4));

                fp=cellfun(@numel,fin_pix_tpj{kk});
                b=find((fp<4));

                ss=cellfun(@mean,pix_tpj{kk},'UniformOutput',false);
                bb=cellfun(@mean,fin_pix_tpj{kk},'UniformOutput',false);

                for si=1:length(s)
                    if ~isempty(s)
                        ss{s(si)}=pix_tpj{kk}{s(si)};
                    end
                end
                for bi=1:length(b)
                    if ~isempty(b)
                        bb{b(bi)}=fin_pix_tpj{kk}{b(bi)};
                    end
                end

                sss=cat(1,ss{:});
                bbb=cat(1,bb{:});
                clear dist_seg Mm Ii
                dist_seg=pdist2(sss,bbb);
                [Mm,Ii] =min(dist_seg,[],1);
                new_t1_rc{kk}=t1_rc{kk}(Ii);

            else
                %new_t1_rc{kk}=[];
                new_t1_rc{kk}=t1_rc{kk};
            end
        end
        act_t1_rc=new_t1_rc;
    
        %%
        if (len_unq_all_tnb==1)
            all_pix_tpj=cat(1,join_pix_tpj{:});
            if (iscell(all_pix_tpj))
                all_pix_tpj=cat(1,all_pix_tpj{:});
            end
            joinpix_tnb= cat(1,joinpix_tnb,pix_tnb{:});
            all_pix_tnb=cat(1,joinpix_tnb{:});
            new_stats_asgn(y_ind).TpPixelList=all_pix_tpj;
            new_stats_asgn(y_ind).TnPixelList=all_pix_tnb;
            new_stats_asgn(y_ind).Tpflag=unq_all_tpj;
            new_stats_asgn(y_ind).Tnflag=unq_all_tnb;
            y_ind= y_ind+1;
            %disp('1st condi');
           
            %*******************************************************************
        elseif (isempty(act_t1_rc) && len_unq_all_tnb>1)
            %disp('halt!');
            witch_sstnb=[];
            sstnb=[];
            newpix_tnb=[];
            joinnewpix_tnb=[];
            sstnb=cell2mat(cellfun(@length,pix_tnb,'uni',0));
            if (numel(sstnb)==1)
                witch_sstnb=find(sstnb==len_unq_all_tnb);
                for w=1:len_unq_all_tnb
                    %put_tnbs=cat(1,put_tnbs,pix_tnb{witch_sstnb}{w});
                    new_stats_asgn(y_ind+w).TpPixelList=fin_pix_tpj{1};
                    new_stats_asgn(y_ind+w).TnPixelList=pix_tnb{witch_sstnb}{w};
                    new_stats_asgn(y_ind+w).Tpflag=unq_all_tpj;
                    %new_stats_asgn(y_ind+w).TnPixelList=pix_tnb{1}{w};
                    new_stats_asgn(y_ind+w).Tnflag=unq_all_tnb(w);
                    %disp('2nd condi part1');
                end

            else
                newpix_tnb = horzcat(pix_tnb{:});
                for s=1:length(unq_all_tnb)
                    indtnb=find(all_tnb==unq_all_tnb(s));
                    joinnewpix_tnb{s}=cat(1,newpix_tnb{indtnb});
                end
                for w=1:len_unq_all_tnb
                    new_stats_asgn(y_ind+w).TpPixelList=fin_pix_tpj{1};
                    new_stats_asgn(y_ind+w).TnPixelList=joinnewpix_tnb{w};
                    new_stats_asgn(y_ind+w).Tpflag=unq_all_tpj;
                    new_stats_asgn(y_ind+w).Tnflag=unq_all_tnb(w);
                    %disp('2nd condi part2');
                end
            end
            y_ind=y_ind+w+1;
            witch_sstnb=[];
            %*******************************************************************
        elseif ((length(act_t1_rc)>1 && len_unq_all_tnb>1) || (length(act_t1_rc)==1 && len_unq_all_tnb>1))
            clear chng_pix_tpj asignt2 chng_centpixtp;
            chng_pix_tpj=fin_pix_tpj;
            chng_centpixtp=centpx_tp;
            done_ind=[];

            for t=1:length(unq_all_tnb)
                t1_rc_tog=[];
                t2_rc_tog=[];
                t1_rc_tog_unq=[];
                tnb_match=[];

                clear tnb_match assign checksame tm finaldist to_lbl;
                clear ap fr ia ib s f;
                ap=1:length(act_t2_rc);
                for na = ap
                    [fr{na},ia{na},ib{na}] = intersect([unq_all_tnb(t)],[act_t2_rc{na}]);
                end
                tnb_match=find(~cellfun('isempty',fr));

                put_t2s=[];
                put_t2scent=[];
                for r=1:length(tnb_match)
                    f=cell2mat(ib(tnb_match(r)));
                    put_t2s=cat(1,put_t2s,pix_tnb{tnb_match(r)}{f});
                end
                put_t2s=sortrows(put_t2s);
                if (numel(put_t2s)==2)
                    put_t2scent=put_t2s;
                else
                    put_t2scent=mean(put_t2s);
                end
                t1_rc_tog=cat(2,act_t1_rc(tnb_match));
                t2_rc_tog=cat(2,act_t2_rc(tnb_match));

                t1_rc_tog_unq=unique(cat(2,t1_rc_tog{:}));
                singleones=[];
                clear len1
                len1=find(cellfun(@length,act_t2_rc)==1);
                singleones= intersect(len1,tnb_match);
                together_lbl={[]};
                to_lbl=[];
                doneids=[];
                non_done=[];

                if (~isempty(singleones))
                    to_lbl=[];
                    together_lbl={};
                    take_t1=[];
                    clear centpx_tmp
                    for tnb_m=1:length(singleones)
                        if  (~isempty(chng_pix_tpj{singleones(tnb_m)}))
                            together_lbl = cat(1,together_lbl,chng_pix_tpj{singleones(tnb_m)});
                        end
                        together_lbl= together_lbl(~cellfun('isempty', together_lbl));
                        chng_pix_tpj{singleones(tnb_m)}={};
                        take_t1=cat(2,take_t1,act_t1_rc{singleones(tnb_m)});
                        act_t1_rc{singleones(tnb_m)}=NaN;
                    end

                    to_lbl=cat(1,together_lbl{:});
                    if (numel(to_lbl)==2)
                        centpx_tmp=to_lbl;
                    else
                        centpx_tmp=mean(to_lbl);
                    end
                    doneids = intersect(singleones,tnb_match);
                    non_done = setxor(tnb_match,doneids);

                    if (~isempty(to_lbl))
                        clear distpq Ip takepix_pj take_t1_tmp;
                        take_t1_tmp=[];
                        takepix_pj={};
                        for tnb_md=1:length(non_done)
                            clear p q
                            p=cell2mat(chng_centpixtp(non_done(tnb_md)));
                            q=centpx_tmp;
                            if  (~isempty(p) && length(p(:,1))>1)
                                [distpq{tnb_md},Ip{tnb_md}]=pdist2(p,q,'euclidean','Smallest',1);
                                takepix_pj{tnb_md}=chng_pix_tpj{non_done(tnb_md)}{Ip{tnb_md}};

                                take_t1_tmp(tnb_md)=act_t1_rc{non_done(tnb_md)}(Ip{tnb_md});
                                chng_pix_tpj{non_done(tnb_md)}{Ip{tnb_md}}=[];
                                chng_centpixtp{non_done(tnb_md)}(Ip{tnb_md},:)=[];

                                for e=1:length(takepix_pj)
                                    if (iscell(takepix_pj{e}))
                                        takepix_pj{e}=cat(1,takepix_pj{e}{:});
                                    end
                                end

                                to_lbl=cat(1,to_lbl,takepix_pj{:});
                                to_lbl= unique(to_lbl,'rows','stable');
                                centpx_tmp=mean(to_lbl);
                                act_t1_rc{non_done(tnb_md)}(Ip{tnb_md})=NaN;

                            elseif (~isempty(p) && ~iscell(p) && length(p(:,1))==1)
                                %f=cell2mat(ib(non_done(tnb_md)));

                                clear newchng_pix_tpj
                                newchng_pix_tpj= chng_pix_tpj{non_done(tnb_md)};
                                takepix_pj{tnb_md}=newchng_pix_tpj;
                                
                                f=find(~cellfun(@isempty,chng_pix_tpj{non_done(tnb_md)}));

                                for s=1:length(takepix_pj)
                                    if (iscell(takepix_pj{s})==0 && isempty(takepix_pj{s})),takepix_pj{s}=[];end
                                    %%if iscell(takepix_pj{s}), takepix_pj{s}=cat(1,takepix_pj{s}{:}); end
                                end
                                chng_pix_tpj{non_done(tnb_md)}={};
                                take_t1_tmp(tnb_md)=act_t1_rc{non_done(tnb_md)}(f);
                                act_t1_rc{non_done(tnb_md)}(f)=NaN;

                                for e=1:length(takepix_pj)
                                    if (iscell(takepix_pj{e}))
                                        takepix_pj{e}=cat(1,takepix_pj{e}{:});
                                    end
                                end
                                to_lbl=cat(1,to_lbl,takepix_pj{:});
                                if iscell(to_lbl),to_lbl=to_lbl{:};end
                                %to_lbl=sortrows(to_lbl);
                                to_lbl= unique(to_lbl,'rows','stable');
                                centpx_tmp=mean(to_lbl);

                            elseif (isempty(p))
                                to_lbl=sortrows(to_lbl);
                                %take_t1_tmp(tnb_md)=[];
                            end
                        end
                        take_t1=unique(cat(2,take_t1,take_t1_tmp));
                    end

                elseif (isempty(singleones) && isempty(to_lbl))
                    clear to_lbltmp
                    for tnb_m=1:length(tnb_match)
                        clear take_t1 assign checksame finaldist ftn gtn tm

                        assign=cell2mat(chng_centpixtp(tnb_match(tnb_m)));
                        checksame=assign-chng_centpixtp{tnb_match(tnb_m)};
                        if (all(checksame==0) |  isnan(checksame(1,1)))
                            if iscell(chng_pix_tpj{tnb_match(tnb_m)})==1
                                ftn=~(cellfun('isempty',chng_pix_tpj{tnb_match(tnb_m)}));
                            else
                                ftn= ~isempty(chng_pix_tpj{tnb_match(tnb_m)});
                            end
                            if sum(ftn)>1
                                [Dt,It]=pdist2(assign,put_t2scent,'euclidean','Smallest',1);
                                to_lbltmp{tnb_m}= chng_pix_tpj{tnb_match(tnb_m)}{It};
                                take_t1=act_t1_rc{tnb_match(tnb_m)}(It);
                                chng_pix_tpj{tnb_match(tnb_m)}{It}=[];
                                act_t1_rc{tnb_match(tnb_m)}(It)=NaN;

                            elseif ((isempty(assign) || isnan(assign(1,1))) && sum(ftn)==0)

                                to_lbltmp{tnb_m}=chng_pix_tpj{tnb_match(tnb_m)};
                                take_t1=[];
                                chng_pix_tpj{tnb_match(tnb_m)}=[];

                            else
                                gtn=find(ftn);
                                Dt=pdist2(assign,put_t2scent,'euclidean','Smallest',1);

                                if (iscell(chng_pix_tpj{tnb_match(tnb_m)})==1)
                                    if ~isempty(gtn)
                                        to_lbltmp{tnb_m}= chng_pix_tpj{tnb_match(tnb_m)}{gtn};
                                        take_t1=act_t1_rc{tnb_match(tnb_m)}(gtn);
                                        chng_pix_tpj{tnb_match(tnb_m)}{gtn}=[];
                                        act_t1_rc{tnb_match(tnb_m)}(gtn)=NaN;
                                    else
                                        to_lbltmp{tnb_m}= chng_pix_tpj{tnb_match(tnb_m)};
                                        take_t1=act_t1_rc{tnb_match(tnb_m)};
                                        chng_pix_tpj{tnb_match(tnb_m)}=[];
                                        %act_t1_rc{tnb_match(tnb_m)}(gtn)=NaN;
                                    end
                                else
                                    to_lbltmp{tnb_m}=chng_pix_tpj{tnb_match(tnb_m)}; %gtn==1
                                    take_t1=act_t1_rc{tnb_match(tnb_m)};
                                    chng_pix_tpj{tnb_match(tnb_m)}=[];
                                    act_t1_rc{tnb_match(tnb_m)}=[];
                                end
                            end
                        else
                            tm =  1:length(tnb_match);
                            for ta =  tm
                                finaldist{ta} = pdist2([assign],[chng_centpixtp{tnb_match(ta)}]);
                            end
                        end
                        if length(tnb_match)==1
                            to_lbl=cat(1,sortrows(to_lbltmp{1}));
                            
                        else
                            to_lbl=(cat(1,(to_lbltmp{:})));
                            if iscell(to_lbl)
                                to_lbl= cat(1,to_lbl{:});
                            end
                        end
                    end
                end

                if exist('doneids')
                    for d=1:length(doneids),act_t2_rc{doneids(d)}=NaN;end % removes only single ids
                end

                for r=1:length(tnb_match)
                    if ~isempty(act_t2_rc{tnb_match(r)})
                        f=cell2mat(ib(tnb_match(r)));
                        act_t2_rc{tnb_match(r)}(f)=NaN;
                    end
                end

                new_stats_asgn(y_ind).TpPixelList=to_lbl;
                new_stats_asgn(y_ind).TnPixelList=put_t2s;
                new_stats_asgn(y_ind).Tpflag=unique(take_t1);
                new_stats_asgn(y_ind).Tnflag=unq_all_tnb(t);
                y_ind=y_ind+1;
            end
        end
        ind_m=[unq_tptn];
        saml_don_id=cat(1,saml_don_id,ind_m);

        %disp('unq_tptn=');disp(unq_tptn);
    end
    new_stats_asgn(all(cell2mat(arrayfun(@(x) structfun(@isempty,x),new_stats_asgn,'UniformOutput',false)),1)) = [];


    clear Assign_tbl A_srttbl
    Assign_tbl = struct2table(new_stats_asgn);
    A_srttbl = sortrows(Assign_tbl,4);
    two_frm_lnks = table2cell(A_srttbl);
    if (~isempty(disappear_t1))
        for z=1:length(disappear_t1)
            append_disapear{z,1}=stats_t1(disappear_t1(z)).PixelList;
            append_disapear{z,2}=[];
            append_disapear{z,3}=disappear_t1(z);
            append_disapear{z,4}=cell2mat([]);
        end
        two_frm_lnks=cat(1,two_frm_lnks,append_disapear);
    end
    Link_list{timeframe}= two_frm_lnks;
end
save(['STIPS_LinkedList','spectmaxr_',num2str(radiuspix),'creatradius_dbscpnts_dbrad_',num2str(creatradius),'_',num2str(dbscanpnts),'_',num2str(dbradius),'cell_centre_75_50_100_100','.mat'],'Link_list','time_inter_sec','px_xysz','Rf1');
