%Aster et al. Nature Communications 2023 Primary Microseism Study Program 2
%displays global trend results
%Requires results from program 1

clear
%list of stations
    station_list=importdata('stations_sel_pub.txt');
    station_list_k={};

%set flags either for acceleration or velocity squared (energy) here
    acc=false;
    vel2=true;

%figure sizing parameters for aspect ratios used in the manuscript
%(three-monitor Mac; may need changing for other computer systems)
 pos=[-2504         128        1784        1130];
 pos1=[-2.5040    0.0460    0.9970    1.2866]*1000;
 pos2=[-2504         128        3700        1130];

    if acc
    unitR='nm s^{-2} y^{-1}';
    unitA='nm s^{-2}';
    titlestr='Acceleration';
    end

    if vel2
    unitR='(nm s^{-1})^2 y^{-1}';
    unitA='(nm s^{-1})^2';
    titlestr='Energy';
    end

%set number of columns in Fig 1/9 plots
    ncols=15;
    req_years=20;
    req_comp=0.8;
    j=1;

%Fig. 1, or Fig. 9, depending on choice of accelration or energy
    figure(300)
    clf

for ista=1:length(station_list)
    clearvars -except coords req_years req_comp...
        j ista station_list station_list_k B_res B1_res B2_res B2h_res B2h2000_res Pct_res Pct1_res ...
        Pct2_res Pct2h_res Pct2h2000_res maxddy minddy acc vel2 unitA unitR titlestr ncols pos* gcol *store

%Primary bands (s)
    band_range='19_17_15';

%load results from first program (generated by Analyze_prim_data_1_Aster_2023.m)
    if acc
    eval(['load Results/results2.0_',band_range,'_',char(station_list(ista)),'_acc.mat']);
    end
    
    if vel2
    eval(['load Results/results2.0_',band_range,'_',char(station_list(ista)),'_vel2.mat']);
    end

    maxddy(j)=max(ddy(~isnan(psd_cull_dy_2)));
    minddy(j)=min(ddy(~isnan(psd_cull_dy_2)));
    ddy_store{j}=ddy(~isnan(psd_cull_dy_2));
    ddy_all_store{j}=ddy;

    psd_cull_dy_2_store{j}=psd_cull_dy_2;
    
    hfun_store{j}=hfun;
    Nyears=(maxddy(j)-minddy(j))/T0;
    Cmp_frac=numel(psd_cull_dy_2(~isnan(psd_cull_dy_2)))/numel(psd_cull_dy_2);

%verify that the required number of years and completeness is satisfied
   if Nyears >= req_years && Cmp_frac >= req_comp

    disp([char(station_list(ista)),' kept: ',num2str(Nyears,3.1),' ',num2str(Cmp_frac),' (',num2str(j),')'])
    B_res(j,:)=[B(1),Bstds(1),B(2),Bstds(2)]*T0;
    B1_res(j,:)=[B1(1),B1stds(1),B1(2),B1stds(2)]*T0;
    B2_res(j,:)=[B2(1),B2stds(1),B2(2),B2stds(2)]*T0;
    B2h_res(j,:)=[B2h(1),B2hstds(1),B2h(2),B2hstds(2)]*T0;
    B2h2000_res(j,:)=[B2h2000(1),B2h2000stds(1),B2h2000(2),B2h2000stds(2)]*T0;
    
    Pct_res(j,:)=[pctperyearB,pctperyearB_std];
    Pct1_res(j,:)=[pctperyearB1,pctperyearB1_std];
    Pct2_res(j,:)=[pctperyearB2,pctperyearB2_std];
    Pct2h_res(j,:)=[pctperyearB2h,pctperyearB2h_std];
    Pct2h2000_res(j,:)=[pctperyearB2h2000,pctperyearB2h2000_std];

    station_list_k(j)=cell(station_list(ista));
    coords(j,:)=str2num(cell2mat(latlon(1)));
    j=j+1;

figure(300)
%plot_smoothing_factor for timeseries (3*T0 for Nature Commun. Paper)
psmooth=round(3*T0);
%psmooth=1;
dsamp=5;
sp_hand1=subplot(ncols,4,j-1);
%     sp_hand1=subplot(ncols,4,i);
    pos1 = get(sp_hand1, 'Position'); % gives the position of current sub-plot
    new_pos1 = pos1 +[0 0 0 0.0005];
    set(sp_hand1, 'Position',new_pos1 ); % set new position of current sub - plot
%exclude seasonal harmonics
smdata=movmedian(psd_cull_dy_2-hfun,psmooth,'omitmissing','Endpoints','shrink');
%include seasonal harmonics
%smdata=movmedian(psd_cull_dy_2,psmooth,'omitmissing','Endpoints','shrink');
dyplot=5*mad(smdata);
ylim([nanmedianB2h-dyplot,nanmedianB2h+dyplot])

smdata(isnan(psd_cull_dy_2-hfun))=NaN;
smdata_store{j-1}=smdata;
ddy10=downsample(ddy,dsamp);
smdata10=downsample(smdata,dsamp);

hold on

ind2000=find(ddy>=datenum(2000,1,1));

xlim([datenum(1987,1,1),datenum(2025,1,1)])
ytickformat('%2.1f')
grid on
grid minor

T(j-1)=title(st_name,'color',[0 0 0],'fontsize',13);
gcol(j-1,:)=[0 0 0];

plot(ddy10,smdata10,'.','markersize',8,'color',[0.5 0.5 0.5])

%define regions for figure station name and plot coloring scheme

%Europe/NE Atlantic Region
if coords(j-1,1)>0 && abs(coords(j-1,2)) < 43
T(j-1)=title(st_name,'color',[0 0 1],'fontsize',13);
plot(ddy,smdata,'.','markersize',1,'color','b')
gcol(j-1,:)=[0 0 1];
end

%Mid-North America Atlantic/Pacific Region
if coords(j-1,1)>0 && coords(j-1,1) < 55 && coords(j-1,2) > -111 && coords(j-1,2) < -60
        T(j-1)=title(st_name,'color',[0.6 0 0],'fontsize',13);
        plot(ddy,smdata,'.','markersize',1,'color',[0.6 0 0])
        gcol(j-1,:)=[0.6 0 0];
end

%Southeast hemisphere
if coords(j-1,1)<0 && coords(j-1,2) > 0
        T(j-1)=title(st_name,'color',[0 0.5 0],'fontsize',13);
        plot(ddy,smdata,'.','markersize',1,'color',[0 0.5 0]')
        gcol(j-1,:)=[0 0.5 0];

end

%Southwest Hemisphere
if coords(j-1,1)<0 && coords(j-1,2) < 0
        T(j-1)=title(st_name,'color',[0 0.7 0.7],'fontsize',13);
        plot(ddy,smdata,'.','markersize',1,'color',[0 0.7 0.7])
    gcol(j-1,:)=[0 0.7 0.7];
end

%declining trends
if B2h_res(j-1,3)+3*B2h_res(j-1,4)<0;set(T(j-1),'fontangle','italic');end


%plot total L1 trend
plot(ddy10,B2h(1)+ddy10*B2h(2),'r','linewidth',4)

%plot post-2000 L1 trend
plot(ddy(ind2000),B2h2000(1)+ddy(ind2000)*B2h2000(2),'k-.','linewidth',3)
xlim([datenum(1987,1,1),datenum(2023,1,1)])

ticks=linspace(floor(min(smdata)),ceil(max(smdata)),5);

set(gca,'ytick',ticks);
    set(gca,'xticklabel',[])
    set(gca,'fontweight','bold')
    set(gca,'fontsize', 10)
    set(T(j-1),'fontsize',14)

else
 disp([char(station_list(ista)),' rejected: ',num2str(Nyears,3),' ',num2str(Cmp_frac)])
hold off
   end
   end

Nsta=j-1;

ticks=datenum(1987:3:2025,1,1);
set(gca,'xtick',ticks);
datetick('x','yyyy','keepticks')
xlabel('Year')
ylabel(unitA)
xlim([datenum(1987,1,1),datenum(2023,1,1)])
set(gca,'fontweight','bold')
set(gca,'fontsize', 10)
set(T(j-1),'fontsize',14)
% 
subplot(ncols,4,Nsta+1)
plot([0 1],[0 1],'color','none')
set(gca,'xtick',[],'ytick',[])
set(gca,'color','none');
set(gca,'XColor','none','YColor','none')
box off
x1=-.2;
y1=1;
dytext=0.5;
%labels at bottom of figure
text(x1,y1,'Eastern N America','color',[0.6 0 0],'fontsize',25)
text(x1,y1-dytext,'N Europe','color','b','fontsize',25)
text(x1,y1-2*dytext,'S and E Hemisphere','color',[0 0.5 0],'fontsize',25)
text(x1,y1-3*dytext,'S and W Hemisphere','color',[0 0.7 0.7],'fontsize',25)
text(x1,y1-4*dytext,'N Pacific','color','k','fontsize',25)
text(x1,y1-5*dytext,'Increasing Trend','color',[0.6 0.6 0.6],'fontsize',25)
text(x1,y1-6*dytext,'Declining Trend','color',[0.6 0.6 0.6],'fontsize',25,'Fontangle','italic')
set (gcf,'position',pos1);

[~,indsort1]=sort(B2h_res(:,3),'descend');

%Blue Marble background for global plot background
load bluemarble.mat

%supplemental figure 
if(acc)
%Global median levels figure (Supp. Fig. 10)
figure(196)
clf
worldmap("world")
geoshow(A,R)
set(gcf,'position',pos)
mfac=17;

hold on
for i=1:Nsta
    med_acc(i)=nanmedian(psd_cull_dy_2_store{i});
    H=plotm(coords(i,1),coords(i,2),'.','color',[0 .8 0],'markersize',mfac*med_acc(i));
    plotm(coords(i,1),coords(i,2),'wo','markersize',mfac*med_acc(i)*.31,'linewidth',2)
    textm(coords(i,1)-.06*mfac*abs(med_acc(i)),coords(i,2),station_list_k(i),'Color','[0 .8 0]','fontsize',20,'fontweight','bold')
end
hold off
title('Primary Microseism Median Acceleration 14 - 20 s','fontsize',40,'fontweight','bold','interpreter','none');
[~,obj]=legend(H,[' 10 ',unitA],'fontsize',35,'location','southoutside');
objhl = findobj(obj, 'type', 'line'); %// objects of legend of type line
set(objhl, 'Markersize', 10*mfac); % set marker size as desired
legend boxoff

set(gcf,'position',pos);
end

%Rate vs. Median Plot (Fig. 4) when acc is selected
    if acc
    figure(198)
    clf
    hold on
    for i=1:52
    nm(i)=nanmedian(psd_cull_dy_2_store{i}-hfun_store{i});
    snm(i)=mad(psd_cull_dy_2_store{i}-hfun_store{i});;
    end
    errorbar(nm,B2h_res(:,3),3*B2h_res(:,4),'o')
    bookfonts
    xlabel('Median Amplitude (nm~s^{-2})','fontsize',24)
    if acc
    ylabel(['R_A (',unitR,')'],'fontsize',24)
    end
    if vel2
    ylabel(['R_E (',unitR,')'],'fontsize',24)
    end
    
    for i=1:52
    text(nm(i),B2h_res(i,3),['  ',station_list_k(i)],'color',gcol(i,:),'fontsize',20)
    plot(nm(i),B2h_res(i,3),'.','color',gcol(i,:),'Markersize',25)
    grid on
    grid minor
    end
    hold on
    plot([0 12],[0 0.005]*12,'k--')
    plot([0 12],-[0 0.005]*12,'k--')
    plot([0 12],[0 0.0025]*12,'k--')
    plot([0 12],-[0 0.0025]*12,'k--')
    plot([0 12],-[0 0]*12,'k--')
    hold off
    end

    %Fig. 3 in manuscript (top or bottom two panels, as selected for acceleration or
    %energy)
    figure(200)
    clf
    subplot(2,1,1)
    ind=1:Nsta;
   
    hold on
    errorbar(ind,B2_res(indsort1(ind),3),3*B2_res(indsort1(ind),4),'.','markersize',15,'linewidth',2,'color','b')
    errorbar(ind,B2h_res(indsort1(ind),3),3*B2h_res(indsort1(ind),4),'.','markersize',20,'linewidth',3,'color','k')
    hold off
    legend('Annual Harmonics Present','Annual Harmonics Removed','fontsize',18,'location','northeast')
    grid on
    ylabel(['R (',unitR,')'])
    yloc=0;
    for i=1:Nsta
        T=text(i,yloc,char(station_list_k(indsort1(i))),'rotation',-55,'fontsize',13);
        if B2h_res(indsort1(i),3)<0
            set(T,'color','[0 0.6 0.6]')
        else
            set(T,'color','r')
        end
    end
    xlim([0,Nsta+1])
   set(gca,'FontSize',20,'LineWidth',1.0);

   subplot(2,1,2)
    [~,indsort2]=sort(Pct2h_res(:,1),'descend');
    ind=1:Nsta;
    hold on
    errorbar(ind,Pct2_res(indsort2(ind),1),3*Pct2_res(indsort2(ind),2),'.','markersize',15,'linewidth',2,'color','b')
    errorbar(ind,Pct2h_res(indsort2(ind),1),3*Pct2h_res(indsort2(ind),2),'.','markersize',20,'linewidth',3,'color','k')
    hold off
    legend('Annual Harmonics Present','Annual Harmonics Removed','fontsize',18,'location','northeast')
    grid on
    ylabel('R')

    yloc=0;

    for i=1:Nsta
        T=text(i,yloc,char(station_list_k(indsort2(i))),'rotation',-55,'fontsize',13);
        if Pct2h_res(indsort2(i),1)<0
            set(T,'color','[0 0.6 0.6]')
        else
            set(T,'color','r')
        end
    end

xlim([0 Nsta+1])
   set(gca,'FontSize',20,'LineWidth',1.0);
xlabel('Trend Rank')
set(gcf,'position',[pos(1) pos(2) pos(3) pos(4)])


%Map Fig. 2 b or c (depending on selection of acceleration or energy
figure(202)
clf
worldmap("world")
geoshow(A,R)
set(gcf,'position',pos)

count3pos=0;
count2pos=0;
count1pos=0;
count0pos=0;
count3neg=0;
count2neg=0;
count1neg=0;
count0neg=0;

mfac=400;
hold on
[~,pctind]=sort(abs(Pct2h_res),'descend');
for i=1:Nsta
    k=pctind(i);

    if Pct2h_res(k,1) > 0

                %three sigma signficance
        if Pct2h_res(k,1)-3*Pct2h_res(k,2) > 0
            H=plotm(coords(k,1),coords(k,2),'r.','markersize',mfac*Pct2h_res(k));
            plotm(coords(k,1),coords(k,2),'wo','markersize',mfac*Pct2h_res(k)*.29,'linewidth',2)
            count3pos=count3pos+1;
        end

                 %two sigma signficance
        if Pct2h_res(k,1)-2*Pct2h_res(k,2) > 0 && Pct2h_res(k,1)-3*Pct2h_res(k,2) < 0
             plotm(coords(k,1),coords(k,2),'ro','markersize',mfac*Pct2h_res(k)*.29,'linewidth',6)
             plotm(coords(k,1),coords(k,2),'wo','markersize',mfac*Pct2h_res(k)*.33,'linewidth',2)
             count2pos=count2pos+1;
        end       

                %one sigma signficance
        if Pct2h_res(k,1)-1*Pct2h_res(k,2) > 0 && Pct2h_res(k,1)-2*Pct2h_res(k,2) < 0
            plotm(coords(k,1),coords(k,2),'ro','markersize',mfac*Pct2h_res(k)*.29,'linewidth',2)
            count1pos=count1pos+1;
        end

        %less than one sigma signficance
        if Pct2h_res(k,1)-1*Pct2h_res(k,2) < 0
            plotm(coords(k,1),coords(k,2),'rs','markersize',5)
            count0pos=count0pos+1;
        end

    textm(coords(k,1)-.06*mfac*abs(Pct2h_res(k)),coords(k,2),station_list_k(k),'Color','r','fontsize',20,'fontweight','bold')
    
    else
        
        %three sigma signficance
        if Pct2h_res(k,1)+3*Pct2h_res(k,2) < 0
            plotm(coords(k,1),coords(k,2),'c.','markersize',abs(mfac*Pct2h_res(k)));
            plotm(coords(k,1),coords(k,2),'wo','markersize',abs(mfac*Pct2h_res(k)*.29),'linewidth',2)
            count3neg=count3neg+1;
        end

         %two sigma signficance
        if Pct2h_res(k,1)+2*Pct2h_res(k,2) < 0 && Pct2h_res(k,1)+3*Pct2h_res(k,2) > 0
             plotm(coords(k,1),coords(k,2),'co','markersize',abs(mfac*Pct2h_res(k)*.29),'linewidth',6)
             plotm(coords(k,1),coords(k,2),'wo','markersize',abs(mfac*Pct2h_res(k)*.33),'linewidth',2)
             count2neg=count2neg+1;
        end  

        %one sigma signficance
        if Pct2h_res(k,1)+1*Pct2h_res(k,2) < 0 && Pct2h_res(k,1)+2*Pct2h_res(k,2) > 0
            plotm(coords(k,1),coords(k,2),'co','markersize',abs(mfac*Pct2h_res(k)*.29),'linewidth',2)
            count1neg=count1neg+1;
        end

        %less than one sigma signficance
        if Pct2h_res(k,1)+1*Pct2h_res(k,2) > 0
            plotm(coords(k,1),coords(k,2),'cs','markersize',5)
            count0neg=count0neg+1;
        end

    textm(coords(k,1)-.06*mfac*abs(Pct2h_res(k)),coords(k,2),station_list_k(k),'Color',[0,0.8,0.8],'fontsize',20,'fontweight','bold')
    end

end
hold off
title(['Microseism Seismic ',titlestr,' Trends ',band_range,' s'],'fontsize',40,'fontweight','bold','interpreter','none');
[~,obj]=legend(H,'  +0.5 %y^{-1} rel. median','fontsize',35,'location','southoutside');
objhl = findobj(obj, 'type', 'line'); %// objects of legend of type line
set(objhl, 'Markersize', 0.5*mfac); % set marker size as desired
legend boxoff

set(gcf,'position',pos);

disp('All data proportional rates:')
disp(['>3-sigma (pos/neg): ',num2str(count3pos),' ',num2str(count3neg)])
disp(['>2-sigma (pos/neg): ',num2str(count2pos),' ',num2str(count2neg)])
disp(['>1-sigma (pos/neg): ',num2str(count1pos),' ',num2str(count1neg)])
disp(['<1-sigma (pos/neg): ',num2str(count0pos),' ',num2str(count0neg)])

figure(205)
clf
worldmap("world")
geoshow(A,R)
set(gcf,'position',pos)

if acc
mfac=600/max(B2h_res(:,3));
legendcirc=round(100*median(abs(B2h_res(:,3))))/100;
end

if vel2
mfac=100;
legendcirc=0.2;
end

hold on
[~,absind]=sort(abs(B2h_res(:,3)),'descend');
for i=1:Nsta
    k=absind(i);

    if B2h_res(k,3) > 0
                %three sigma signficance
        if B2h_res(k,3)-3*B2h_res(k,4) > 0
            H=plotm(coords(k,1),coords(k,2),'r.','markersize',mfac*B2h_res(k,3));
            plotm(coords(k,1),coords(k,2),'wo','markersize',mfac*B2h_res(k,3)*.29,'linewidth',2)
        end

                 %two sigma signficance
        if B2h_res(k,3)-2*B2h_res(k,4) > 0 && B2h_res(k,3)-3*B2h_res(k,4) < 0
             plotm(coords(k,1),coords(k,2),'ro','markersize',mfac*B2h_res(k,3)*.29,'linewidth',6)
             plotm(coords(k,1),coords(k,2),'wo','markersize',mfac*B2h_res(k,3)*.33,'linewidth',2)
        end       

                %one sigma signficance
        if B2h_res(k,3)-1*B2h_res(k,4) > 0 && B2h_res(k,3)-2*B2h_res(k,4) < 0
        plotm(coords(k,1),coords(k,2),'ro','markersize',mfac*B2h_res(k,3)*.29,'linewidth',2)
        end

                %less than one sigma signficance
        if B2h_res(k,3)-1*B2h_res(k,4) < 0
        plotm(coords(k,1),coords(k,2),'rs','markersize',5)
        end

    %textm(coords(k,1)-.06*mfac*abs(B2h_res(k,3)),coords(k,2),station_list_k(k),'Color','k','fontsize',20)
    textm(coords(k,1)-.06*mfac*abs(B2h_res(k,3)),coords(k,2),station_list_k(k),'Color','r','fontsize',20,'fontweight','bold');
    
    else
        
        %three sigma signficance
        if B2h_res(k,3)+3*B2h_res(k,4) < 0
            plotm(coords(k,1),coords(k,2),'c.','markersize',abs(mfac*B2h_res(k,3)));
            plotm(coords(k,1),coords(k,2),'wo','markersize',abs(mfac*B2h_res(k,3)*.29),'linewidth',2)
        end

         %two sigma signficance
        if B2h_res(k,3)+2*B2h_res(k,4) < 0 && B2h_res(k,3)+3*B2h_res(k,4) > 0
             plotm(coords(k,1),coords(k,2),'co','markersize',abs(mfac*B2h_res(k,3)*.29),'linewidth',6)
             plotm(coords(k,1),coords(k,2),'wo','markersize',abs(mfac*B2h_res(k,3)*.33),'linewidth',2)
        end  

        %one sigma signficance
        if B2h_res(k,3)+1*B2h_res(k,2) < 0 && B2h_res(k,3)+2*B2h_res(k,4) > 0
        plotm(coords(k,1),coords(k,2),'co','markersize',abs(mfac*B2h_res(k,3)*.29),'linewidth',2)
        end

        %less than one sigma signficance
        if B2h_res(k,3)+1*B2h_res(k,4) > 0
        plotm(coords(k,1),coords(k,2),'cs','markersize',5)
        end

    textm(coords(k,1)-.06*mfac*abs(B2h_res(k,3)),coords(k,2),station_list_k(k),'Color',[0,0.8,0.8],'fontsize',20,'fontweight','bold');
    end

end
hold off
title(['Microseism Seismic ',titlestr,' Trends ',band_range,' s'],'fontsize',40,'fontweight','bold','interpreter','none');
[~,obj]=legend(H,['  +',num2str(legendcirc),' ',unitR],'fontsize',35,'location','southoutside');
objhl = findobj(obj, 'type', 'line'); %// objects of legend of type line
set(objhl, 'Markersize', legendcirc*mfac); % set marker size as desired
legend boxoff

set(gcf,'position',pos);

%Operational Station Histories plotl (Supp. Fig. 9)
figure(203)
clf
hold on
for j=1:Nsta
meandate(j)=mean([maxddy(j),minddy(j)]);
negdate(j)=meandate(j)-minddy(j);
posdate(j)=maxddy(j)-meandate(j);
%reverse index to produce alphabetical listing from the top down
k=Nsta-j+1;
ddy=cell2mat(ddy_store(j));
ddy10=downsample(ddy,10);
optime10=k*ones(size(ddy10));
lineH=plot(ddy10,optime10,'.','markersize',5);
color = get(lineH, 'Color');
text(datenum(1985,10,1),k,char(station_list_k(j)),'fontsize',18,'color',color,'fontweight','bold')
end
hold off
grid on
ylim([0,Nsta+1])
bookfonts
datetick('x',10)
set(gca,'ytick',[])
xlim([datenum(1985,1,1),datenum(2024,1,1)])
set(gcf,'position',[pos(1) pos(2) pos(3)/2 pos(4)])
bookfonts
grid on

%Supp. Fig. 12 or 13 depending of choice of acceleration or energy
figure(304)
clf
for i=1:Nsta
    hold on
errorbarxy(Pct2h_res(i,1),Pct2h2000_res(i,1),3*Pct2h_res(i,2),3*Pct2h2000_res(i,2),{'ko','k','k'});
plot(Pct2h_res(i,1),Pct2h2000_res(i,1),'.','color',gcol(i,:),'markersize',30)
    T=text(Pct2h_res(i,1),Pct2h2000_res(i,1),['   ',char(station_list_k(i))],'fontsize',15,'fontweight','bold','rotation',-45);
        set(T,'color',gcol(i,:))
end
plot([-0.6 1.4],[-0.6 1.4],'-.')
hold off
if acc
xlabel('P_A (% y^{-1}; All)')
ylabel('P_A (% y^{-1}; Post-2000)')
end
if vel2
xlabel('P_E (% y^{-1}; All)')
ylabel('P_E (% y^{-1}; Post-2000)')   
end
grid on
ylim([-.6 1.4]);
xlim([-.6 1.4]);
axis square
bookfonts

disp(['Global mean trends (abs, pct)y^{-1}: ',num2str(mean(B2h_res(:,3))),' ',num2str(mean(Pct2h_res(:,1)))])
disp(['Global Average 3 sigmas: ',num2str(sqrt(mean(B2h_res(:,4)).^2)),', ',num2str(sqrt(mean(Pct2h_res(:,2)).^2))])

%save workspace for use in Program 3
if acc
eval(['save results_workspace_acc_',band_range,'.mat *store coords Pct2* B2* station_list_k gcol deltad Nsta']);
end 
if vel2
eval(['save results_workspace_vel2_',band_range,'.mat *store coords Pct2* B2* station_list_k gcol deltad Nsta']);
end





