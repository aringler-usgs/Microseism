%Aster et al. Nature Communications 2023 Primary Microseism Study Program 3
%
%Requres m_map mapmaking software
% (Pawlowicz, R., 2020. "M_Map: A mapping package for MATLAB", version 1.4m; www.eoas.ubc.ca/~rich/map.html)
% and distinguishable_colors
%(Tom Holy; https://www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors)
%(included in this distribution)
%
%Requires results from program 2

clear
acc=true;
vel2=false;
unitA='nm s^{-2}';
unitE='Energy';

    %choosing acc generates cluster and other acc-derived plots
    if acc
    load results_workspace_acc_19_17_15.mat
    units=unitA;
    end

    %choosing vel2 generates the velocity squared global stack plot and correlation
    %for various smoothing lengths
    if vel2
    units=unitE;
    load results_workspace_vel2_19_17_15.mat
    end

T0=365.242;
Nsta=length(ddy_all_store);

Nsta=length(station_list_k);
for i=1:Nsta
    ddy=cell2mat(ddy_all_store(i));
    smdata=(psd_cull_dy_2_store{i}-hfun_store{i});
    %remove trends for correlation clustering
    smdata_dt=smdata-(B2h_res(i,1)/T0+ddy*B2h_res(i,3)/T0);
    smdata_dt_sm1=movmedian(smdata_dt,T0/6,'omitmissing','endpoints','fill');
    smdata_dt_store{i}=smdata_dt;
    smdata_dt_store_sm1{i}=smdata_dt_sm1;
end

%introduce climate indices for clustering
n_ind=2;
load Indices/SOI.mat
load Indices/BEST.mat
%load Indices/PDO.mat
%load Indices/WP.mat
%load Indices/NAO.mat
%increase "Station" indices to include climate indices in correlations here
Nsta=Nsta+n_ind;

smdata_dt_store_sm1{end+1}=detrend(SOI_data_interp);
smdata_dt_store_sm1{end+1}=detrend(BEST_data_interp);
%smdata_dt_store_sm1{end+1}=detrend(NAO_data_interp);
%smdata_dt_store_sm1{end+1}=detrend(PDO_data_interp);
%smdata_dt_store_sm1{end+1}=detrend(WP_data_interp);
%smdata_dt_store_sm1{end+1}=detrend(NAO_data_interp);
ddy_all_store{end+1}=dinterp;
ddy_all_store{end+1}=dinterp;
ddy_all_store{end+1}=dinterp;
ddy_all_store{end+1}=dinterp;
ddy_all_store{end+1}=dinterp;
ddy_all_store{end+1}=dinterp;
station_list_k{end+1}='\bf SOI';
station_list_k{end+1}='\bf ENSO';
%station_list_k{end+1}='\bf NAO';
%station_list_k{end+1}='\bf PDO';
%list_k{end+1}='\bf WPI';
%station_list_k{end+1}='\bf AMO';
for i=1:Nsta
    for j=1:i
        x1=cell2mat(smdata_dt_store_sm1(i));
        x2=cell2mat(smdata_dt_store_sm1(j));
    x1(isnan(x1))=0;
    x2(isnan(x2))=0;
    x1=x1-mean(x1);
    x2=x2-mean(x2);
    norm=sqrt(sum(x1.^2)*sum(x2.^2));
    xc=xcorr(x1,x2)/norm;

    %correct for time index lag between two stations
    [cr,lr]=xcorr(ddy_all_store{i},ddy_all_store{j});
    [cmax,icmax]=max(cr);
   C(i,j)=xc(icmax);
    end
end
C=(C+tril(C,-1)');

figure(23)
%Correlation matrix figure (Supp. Fig. 14)
imagesc(C)
axis xy
colorbar
set(gca,'ytick',1:Nsta);
set(gca,'xtick',1:Nsta);
set(gca,'yticklabel',station_list_k)
set(gca,'xticklabel',station_list_k)
title(['Primary Microseism (',units,')'])
bookfonts
axis square
bookfonts

%create the "dissimilarity" matrix;
D=1-C;
Z=linkage(D,'ward');

if acc
%dendrogram figure (bottom of Fig. 5)
figure(104)

[H,T,outperm]=dendrogram(Z,Nsta,'colorthreshold',3,'orientation','right');
set(H,'LineWidth',2);
xlabel('D')
set(gca,'yticklabel',station_list_k(outperm))
xlim([0 4])

    title('Primary Microseism')

bookfonts

a=findobj(H,'type','axe');
x=get(get(a,'xlabel'),'string');
end

Nsta=Nsta-n_ind;

if acc
ticklabels_new=station_list_k;
for i=1:Nsta
    c_cell{i}='{0 0 0 }';

%Generic global regions
%Europe/NE Atlantic Region
if coords(i,1)>0 && abs(coords(i,2)) < 43; c_cell{i}='{0 0 1} '; end
%NA/NW Atlantic Region
if coords(i,1)>0 && coords(i,1) < 55 && coords(i,2) > -111 && coords(i,2) < -60; c_cell{i}='{0.6 0 0} '; end
%SE Hemisphere
if coords(i,1)<0 && coords(i,2) > 0; c_cell{i}='{0 0.5 0} '; end
%SW Hemisphere
if coords(i,1)<0 && coords(i,2) < 0; c_cell{i}='{0 0.7 0.7} '; end

    ticklabels_new{i}=['\color[rgb]',c_cell{i},ticklabels_new{i}];
end
set(gca,'yticklabel',ticklabels_new(outperm));

for i=1:Nsta
    for j=1:Nsta
    [azi(i,j),bazi(i,j),range(i,j),delta(i,j)]=edist(coords(i,1),coords(i,2),coords(j,1),coords(j,2));
    end
end

vdel=delta(tril(true(size(delta))));
CC=C(1:52,1:52);
vc=CC(tril(true(size(CC))));

ind=find(vdel>.1);
vc=vc(ind);
vdel=vdel(ind);

%correlation vs distance figure
%Reproduced in Fig. 
figure(105)
vdel_p=vdel;
vc_p=vc;
plot(vdel_p,vc_p,'.','color',[0.8 0 0]','markersize',15)
xlabel('\Delta (degrees)')
ylabel('Correlation')

deld=15;
d=0:deld/2:180-deld;
for i=1:length(d)
    ind=find(vdel_p>=d(i) & vdel_p<d(i)+deld);
    mC_p(i)=mean(vc_p(ind));
    stdC_p(i)=std(vc_p(ind));
end
hold on
plot(d+deld/2,mC_p,'k-','linewidth',5)
plot(d+deld/2,mC_p+stdC_p,'k--','linewidth',5)
plot(d+deld/2,mC_p-stdC_p,'k--','linewidth',5)
hold off
grid on
grid minor
bookfonts

end

for i=1:Nsta
    ddy_min(i)=min(ddy_all_store{i});
    ddy_max(i)=max(ddy_all_store{i});
end
ddy_base=(min(ddy_min):max(ddy_max))';
for i=1:Nsta
    psd2{i}=(psd_cull_dy_2_store{i}-hfun_store{i});
    psd2{i}=psd2{i}/nanmedian(psd2{i});
    data_base_dharm_prop(:,i)=interp1(ddy_all_store{i},psd2{i},ddy_base,'linear',NaN);
end
global_med_prop=nanmedian(data_base_dharm_prop,2);

%maximally distinguishable colors
colors=distinguishable_colors(Nsta+n_ind);

%global stack figure with various (1/6- 3- and 5-year smoothing (3-year smoothing for acceleration and velocity are
%reproduced in Fig. 7) and 1/6-year smoothing is reproduced in Fig. 6.
figure(107)
clf
%Reguero et al., 2009 wave energy increase for plotting
years=1945:2023;
indy=find(years>=1988 & years <= 2023);
years2=years(indy);
dp2=datenum(years2,1,1);

p(1)=2.2;
for i=2:length(years)
    p(i)=p(i-1)*(1+0.47/100);
end
p2=p(indy);
p2=p2/median(p2);

%normalize trends to mid-interval for plotting clarity
p2=p2-(p2(find(years2==2005)))+1-.015;

%All data trend from this study for plotting
pa(1)=1;
for i=2:length(years)
    pa(i)=pa(i-1)+pa(i-1)*(0.27/100);
end
p2a=pa(indy);
p2a=p2a/median(p2a);
%normalize to mid-interval
p2a=p2a-(p2a(find(years2==2005)))+1;

%post-2000 data trend from this study for plotting
pa2000(1)=1;
for i=2:length(years)
pa2000(i)=pa2000(i-1)+pa2000(i-1)*(0.35/100);
end
p2a2000=pa2000(indy);
p2a2000=p2a2000/median(p2a2000);
%normalize to mid-interval
p2a2000=p2a2000-(p2a2000(find(years2==2005)))+1+.015;

smyear=[1/6,3 5 10];
yl=[.4 1.6; .75 1.25; .85 1.15; .92 1.08];

for i=1:3
    subplot(3,1,i)
sm=round(T0*smyear(i));
set(gca,'ColorOrderIndex',1)

all_meds(:,1:Nsta)=movmedian(data_base_dharm_prop(:,1:Nsta),sm,'omitmissing');

for j=1:Nsta
first_data_k(j)=find(~isnan(data_base_dharm_prop(:,j)),1);
all_meds(1:first_data_k(j)+round(sm/2),j)=NaN;
all_meds(end-round(sm/2):end,j)=NaN;
end

hold on
for kk=1:Nsta
plot(ddy_base,all_meds(:,kk),'linewidth',3,'color',[colors(kk,:),0.4])
end
if i==1
    title(['Primary Microseism (',units,')'])
end

hold on
plot([ddy_base(1),ddy_base(end)],[1,1],'-.','color',[0.3 0.3 0.3],'linewidth',3)

%moving median of the median
z=movmedian(nanmedian(data_base_dharm_prop')',sm,'omitmissing');

z(1:min(first_data_k)+round(sm/2))=NaN;
z(end-round(sm/2):end)=NaN;
plot(ddy_base,z,'-','linewidth',5,'color',[0 0 0 0.8])

ind=find(dp2 < datenum(2009,1,1));
plot(dp2(ind),p2(ind),'r-.','linewidth',6)

ind=find(dp2>=datenum(2000,1,1) & dp2 < datenum(2022,10,1));
plot(dp2(ind),p2a2000(ind),'b-.','linewidth',6)

ind=find(dp2>=datenum(1989,1,1) & dp2 < datenum(2022,10,1));
plot(dp2(ind),p2a(ind),'g-.','linewidth',6)

grid on
grid minor
datetick('x')
ylabel([units,' rel. med. (S=',num2str(smyear(i),2),' yr)'])
bookfonts
ylim(yl(i,:))
xlim([datenum(1988,1,1) datenum(2023,1,1)])
ylrange=yl(i,2)-yl(i,1);
end
hold on

hold off

axP = get(gca,'Position'); 
LH=legend([station_list_k(1:Nsta)],...
    'location','southoutside','orientation','horizontal','fontsize',12);
LH.NumColumns=11;
set(gca, 'Position', axP);
set(gcf,'position',[1219 34 1084 1281]);

if acc
figure(115)
clf

smyear=[1/6,3 5 10];
yl=[.4 1.6; .75 1.25; .85 1.15; .92 1.08];

for i=1:3
    subplot(3,1,i)
sm=round(T0*smyear(i));
set(gca,'ColorOrderIndex',1)
all_meds(:,1:Nsta)=movmedian(data_base_dharm_prop(:,1:Nsta),sm,'omitmissing');

for j=1:Nsta
first_data_k(j)=find(~isnan(data_base_dharm_prop(:,j)),1);
all_meds(1:first_data_k(j)+round(sm/2),j)=NaN;
all_meds(end-round(sm/2):end,j)=NaN;
end

hold on
for kk=1:Nsta
plot(ddy_base,all_meds(:,kk),'linewidth',3,'color',[colors(kk,:),0.4])
end
if i==1
    title('Primary Microseism')
end
hold on
plot([ddy_base(1),ddy_base(end)],[1,1],'-.','color',[0.3 0.3 0.3],'linewidth',3)

%moving median of the median
z=movmedian(nanmedian(data_base_dharm_prop')',sm,'omitmissing');

z(1:min(first_data_k)+round(sm/2))=NaN;
z(end-round(sm/2):end)=NaN;
plot(ddy_base,z,'-','linewidth',5,'color',[0 0 0 0.8])


grid on
grid minor
datetick('x')
ylabel([units,' rel. med. (S=',num2str(smyear(i),2),' yr)'])
bookfonts
ylim(yl(i,:))
xlim([datenum(1988,1,1) datenum(2023,1,1)])
ylrange=yl(i,2)-yl(i,1);

end

axP = get(gca,'Position'); 
set(gca, 'Position', axP);
set(gcf,'position',[1 81 476 1281]);

%regional stacks across station clusters
    D2=1-C;
    Z2=linkage(D2,'ward');
ii=cluster(Z2,7);

%1/6-year smoothing (reproduced in Fig. 6)
figure(500)
clf

%3-year smoothing (reproduced in Fig. 7)
figure(501)
clf

%loop over the seven clusters
for l=1:7
%find indices for this cluster
indrc{l}=find(ii==l);
indr=indrc{l};

%produce colleciotns of each cluster with various degrees (1/6, 3, 5-y) of
%smoothing
figure(107+l)
clf
%loop over various smoothing windows
for i=1:3
    subplot(3,1,i)
sm=round(T0*smyear(i));

%moving medians of each time series
reg_meds=movmedian(data_base_dharm_prop(:,indr(indr<53)),sm,'omitmissing');
for j=1:length(indr(indr<53))
    first_data_i(j)=find(~isnan(data_base_dharm_prop(:,indr(j))),1);
reg_meds(1:first_data_i(j)+round(sm/2),j)=NaN;
reg_meds(end-round(sm/2):end,j)=NaN;
end

%moving median of the median
z=movmedian(nanmedian(data_base_dharm_prop(:,indr(indr<53))')',sm,'omitmissing');
    z(1:round(sm/2))=NaN;
    z(end-round(sm/2):end)=NaN;
    z(1:min(first_data_i)+round(sm/2))=NaN;

hold on
for kk=1:length(indr(indr<53))
plot(ddy_base,reg_meds(:,kk),'linewidth',3,'color',[colors(indr(kk),:),0.5])
end

    if indr(end)==53
    soi=movmedian(SOI_data_interp,sm);
    soi=soi-mean(soi);
    soi=(soi/max(soi))*3*mad(z)+1;
    soi_base=interp1(dinterp,soi,ddy_base);
    soi_base(1:round(sm/2)+min(first_data_i))=NaN;
    soi_base(end-round(sm/2):end)=NaN;
    plot(ddy_base,soi_base,'r-','linewidth',5)
    end

    if indr(end)==54
    best=movmedian(BEST_data_interp,sm);
    best=best-mean(best);
    best=(best/max(best))*3*mad(z)+1;
    best_base=interp1(dinterp,best,ddy_base);
    best_base(1:round(sm/2)+min(first_data_i))=NaN;
    best_base(end-round(sm/2):end)=NaN;
    plot(ddy_base,best_base,'r-','linewidth',5)
    end

if i==1
    title(['Primary Microseism (Cluster ',num2str(l),'p)'])
end

hold on
plot([ddy_base(1),ddy_base(end)],[1,1],'-.','color',[0.3 0.3 0.3],'linewidth',3)

plot(ddy_base,z,'-','linewidth',4,'color',[0 0 0 0.8])

grid on
grid minor
datetick('x')
ylabel([units,' rel. med. (S=',num2str(smyear(i),2),' yr)'])
bookfonts

ylim(yl(i,:))

xlim([datenum(1988,1,1) datenum(2023,1,1)])
ylrange=yl(i,2)-yl(i,1);
if i==1
    figure(500)
    set(gcf,'Position',[-2017 19 1704 1343])
    subplot(4,2,l)
    hold on

    for kk=1:length(indr(indr<53))
    plot(ddy_base,reg_meds(:,kk),'linewidth',3,'color',[colors(indr(kk),:),0.5])
    end

    if indr(end)==53
    soi=movmedian(SOI_data_interp,sm);
    soi=soi-mean(soi);
    soi=(soi/max(soi))*3*mad(z)+1;
    soi_base=interp1(dinterp,soi,ddy_base);
    soi_base(1:round(sm/2)+min(first_data_i))=NaN;
    soi_base(end-round(sm/2):end)=NaN;
    plot(ddy_base,soi_base,'r-','linewidth',5)
    end

    if indr(end)==54
    best=movmedian(BEST_data_interp,sm);
    best=best-mean(best);
    best=(best/max(best))*3*mad(z)+1;
    best_base=interp1(dinterp,best,ddy_base);
    best_base(1:round(sm/2)+min(first_data_i))=NaN;
    best_base(end-round(sm/2):end)=NaN;
    plot(ddy_base,best_base,'r-','linewidth',5)
    end

    H=gca;
    LH=legend([station_list_k(indr),'',''],'location','eastoutside','orientation','horizontal',...
    'fontsize',13);
    LH.NumColumns=1;
    ylim([yl(i,1),yl(i,2)])
    datetick('x')
    
    ylabel([units,' rel. med.'])
  
    title(['Cluster ',char(station_list_k(indr(1))),' (S=',num2str(smyear(i),2),' yr)'],'fontsize',20);

     bookfonts
     grid on
     grid minor
     hold on
     plot([ddy_base(1),ddy_base(end)],[1,1],'-.','color',[0.3 0.3 0.3],'linewidth',3,'HandleVisibility','off')
     plot(ddy_base,z,'k-','linewidth',4,'HandleVisibility','off','color',[0 0 0 0.8])
     hold off
     set(gca, 'Position', H.Position);
end

if i==3
    figure(501)
    set(gcf,'Position',[-2017 19 1704 1343])
    subplot(4,2,l)
    hold on

    for kk=1:length(indr(indr<53))
    plot(ddy_base,reg_meds(:,kk),'linewidth',3,'color',[colors(indr(kk),:),0.5])
    end

        if indr(end)==53
    soi=movmedian(SOI_data_interp,sm);
    soi=soi-mean(soi);
    soi=(soi/max(soi))*3*mad(z)+1;
    soi_base=interp1(dinterp,soi,ddy_base);
    soi_base(1:round(sm/2)+min(first_data_i))=NaN;
    soi_base(end-round(sm/2):end)=NaN;
    plot(ddy_base,soi_base,'r-','linewidth',5)
    end

    if indr(end)==54
    best=movmedian(BEST_data_interp,sm);
    best=best-mean(best);
    best=(best/max(best))*3*mad(z)+1;
    best_base=interp1(dinterp,best,ddy_base);
    best_base(1:round(sm/2)+min(first_data_i))=NaN;
    best_base(end-round(sm/2):end)=NaN;
    plot(ddy_base,best_base,'r-','linewidth',5)
    end

    H=gca;
    LH=legend([station_list_k(indr),'',''],'location','eastoutside','orientation','horizontal',...
    'fontsize',13);
    LH.NumColumns=1;
    ylim([yl(i,1),yl(i,2)])
    datetick('x')
    
    ylabel([units,' rel. med.'])
   
    title(['Cluster ',char(station_list_k(indr(1))),' (S=',num2str(smyear(i),2),' yr)'],'fontsize',20);

     bookfonts
     grid on
     grid minor
     hold on
     plot([ddy_base(1),ddy_base(end)],[1,1],'-.','color',[0.3 0.3 0.3],'linewidth',3,'HandleVisibility','off')
     plot(ddy_base,z,'-','linewidth',4,'HandleVisibility','off','color',[0 0 0 0.8])
     hold off
     set(gca, 'Position', H.Position);
end
figure(107+l)

end

hold off

axP = get(gca,'Position'); 
LH=legend(station_list_k(indr),'location','southoutside','orientation','horizontal',...
    'fontsize',14);
LH.NumColumns=5;
set(gca, 'Position', axP);
set(gcf,'position',[1 81 476 1281]);

end

%set range of ocean topography to color
dep_range=-300;
%Map figure showing ocean floor and cluster coloring to match dendrogram
%(top of Fig. 5a)
figure(1000)
clf

H=m_proj('robinson','lon',[-180 180]);

set(gcf,'color','w')   % Set background colour before m_image call

c_mat = [74 80 255; 255 0 0; 0 145 255; 73 255 0; 255 0 219; 0 255 145; 255 219 0]/255;

clim([dep_range 0]);
c=colormap(flipud([flipud(m_colmap('blues',10));m_colmap('jet',118)]));
c(1,:)=[0.9 0.9 0.9]*.3;
colormap(c)
m_etopo2('shadedrelief','gradient',3);
 
m_gshhs('c','patch',[.95 .95 .95]);
 
m_grid('box','fancy','fontsize',20);

%show general geographic divisions on map
m_line([-180;180],[0;0],'linewidth',5,'color','k');
m_line(-[43,43].*ones(90,2),(0:89)','linewidth',5,'color','b');
m_line([43,43].*ones(90,2),(0:89)','linewidth',5,'color','b');
m_line([-43;43],[1;1],'linewidth',5,'color','b');
m_line(-[111,111].*ones(56,2),(0:55)','linewidth',5,'color',[143 20 2]/255);
m_line(-[60,60].*ones(56,2),(0:55)','linewidth',5,'color',[143 20 2]/255);
m_line(-[111;60],[55;55],'linewidth',5,'color',[143 20 2]/255);
m_line(-[111;60],[1;1],'linewidth',5,'color',[143 20 2]/255);
m_line(-[180;0],[0;0],'linewidth',5,'color','c');
m_line([180;0],[0;0],'linewidth',5,'color',[2 171 46]/255);
m_line(-[179.5,179.5].*ones(90,2),(-89:0)','linewidth',5,'color','c');
m_line(-[0.75,0.75].*ones(90,2),(-89:0)','linewidth',5,'color','c');
m_line([179.5,179.5].*ones(90,2),(-89:0)','linewidth',5,'color',[2 171 46]/255);
m_line([0.75,0.75].*ones(90,2),(-89:0)','linewidth',5,'color',[2 171 46]/255);

for ii=1:7
          indr=indrc{ii};
          indr=indr(indr<53);
          m_line(coords(indr,2),coords(indr,1),'marker','o','color','w','linewi',2,...
          'linest','none','markersize',25,'markerfacecolor',c_mat(ii,:));
          m_text(coords(indr,2),coords(indr,1),station_list_k(indr),'color',c_mat(ii,:),...
              'FontSize',25)
end

ax=m_contfbar(1.03,[.2 0.8],[dep_range 0],dep_range:5:0,'edgecolor','none','endpiece','no');
xlabel(ax,'Depth (m)','color','k');

set(gcf,'position',[270 545 1599 725]);
set(ax,'FontSize',25)

%ENSO/SOI plots
%Supp. Fig. 15
figure(21)
subplot(2,1,1)
plot(dinterp,SOI_data_interp,'linewidth',2)
hold on
plot(dinterp,BEST_data_interp,'linewidth',2)

datetick
legend('SOI','ENSO')
bookfonts
ylabel('Index')
hold off
xlim([datenum(1988,1,1) datenum(2023,1,1)])

subplot(2,1,2)
sm=3*T0;
    soi=movmedian(SOI_data_interp,sm);
    soi(1:sm)=NaN;
    soi(end-sm:end)=NaN;

    soi_basesm=interp1(dinterp,soi,ddy_base);
    
    best=movmedian(BEST_data_interp,sm);
    best(1:sm)=NaN;
    best(end-sm:end)=NaN;

    best_basesm=interp1(dinterp,best,ddy_base);
plot(ddy_base,soi_basesm,'linewidth',2)
hold on
plot(ddy_base,best_basesm,'linewidth',2)
datetick
legend('SOI','ENSO')
bookfonts
ylabel('Index')
hold off
xlim([datenum(1988,1,1) datenum(2023,1,1)])
end
