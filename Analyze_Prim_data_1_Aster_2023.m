%Aster 2023 Primary Microseism Code Program 1
%Calculates relevant time series and trends from spectral data

clear

%choose period bands to sum
%SF Primary
usebands=1:3;

%tropical year in days for harmonic function estimation
T0=365.242;

%toggle to show diagnostic figures for each station as the code runs
showfigs=true;

%toggle here to calculate and process either acceleration or velocity-squared
%data (choose one of the other)
acc=true;
vel2=false;

if acc
    unitR='nm s^{-2}y^{-1}';
    unitA='nm s^{-2}';
end

if vel2
    unitR='(nm s^{-1})^2y{-1}';
    unitA='(nm s^{-1})^2';
end

%Mid-band periods for each 2 s-width integration band in data file
psd_periods=[19,17,15,13,11,9,7,5];

%Run all stations using this list
stations=importdata('stations_sel_pub.txt');
%Individual or selected station run
stations={'CHTO'};

%get master station information from this file
xstation_master=importdata('stations_total_all.txt');
for i=1:length(xstation_master)
stns_master(i,:)=strsplit(char(xstation_master(i)),',');
end

%load COMCAT earthquake catalogue for large earthquake culling
X_eq=importdata('Earthquake_Global_PDE/M_gt_5.75_short.txt',',',1);

%process each station sequentialy within this loop
for ista=1:length(stations)
clearvars -except ista stations stns_master usebands psd_periods showfigs T0 acc vel2 unitA unitR X_eq
st_name=char(stations(ista));
ind=find(strcmp(stns_master(:,2),st_name) & ~contains(stns_master(:,1),'X'),1,'last');
NET=char(stns_master(ind,1));
GLOC=char(stns_master(ind,6));

%import the PSD station file
%load data (spectral acceleration PSD structure: Z, and associated sorted dates, dates_sort)
load(['PSD_data_mat/',st_name,'_spectra.mat'])

%retrieve data in nm/s^2 (square root of the integrated PSD) into psd_data
%these indices contain the PSD estimates used in this study
psd_data=Z.data(:,[5,7,9,11,13,15,17,19]);

%end

if vel2
%convert data to velocity squared spectra if vel2 is true
for i=1:8
    psd_data(:,i)=(psd_data(:,i)*(psd_periods(i)/(2*pi)));
end
end

%Duplicate raw data for outlier and earthquake culling
psd_raw=zeros(size(psd_data(:,1)));

for iband=usebands
psd_raw=psd_raw+psd_data(:,iband);
end
bstr=sprintf('%.0f ',psd_periods(usebands));
bstrf=sprintf('%.0f_',psd_periods(usebands));

psd=psd_raw;

latlon=Z.textdata(1,1);

%Extract relevant earthquake data from COMCAT Catalogue
Xc=char(X_eq.textdata(2:end,1));
eq_yr=str2num(Xc(:,1:4));
eq_mon=str2num(Xc(:,6:7));
eq_day=str2num(Xc(:,9:10));
eq_hour=str2num(Xc(:,12:13));
eq_min=str2num(Xc(:,15:16));
eq_sec=str2num(Xc(:,18:23));
eq_dnum=datenum([eq_yr,eq_mon,eq_day,eq_hour,eq_min,eq_sec]);

%indices of earthquakes in each of the below COMCAT magnitude ranges
ind0=find(X_eq.data(:,4)>=5.75 & X_eq.data(:,4) < 7);
ind1=find(X_eq.data(:,4)>=7 & X_eq.data(:,4) < 7.5);
ind2=find(X_eq.data(:,4)>=7.5 & X_eq.data(:,4) < 8);
ind3=find(X_eq.data(:,4)>=8 & X_eq.data(:,4) < 9);
ind4=find(X_eq.data(:,4)>=9);

%Large earthquake window rejection
wlen=3.0; %culling window parmater in hours

%1 x wlen culling for this magnitude range
N6=length(ind0);
    ds=eq_dnum(ind0);
for i=1:N6
    [~,ind0_store(i,1)]=min(abs(dates_sort-ds(i)));
end
ind0_eq_w=ind0_store(ind0_store>1);

i_all0=[];

for j=1:length(ind0_eq_w)
    irange=ind0_eq_w(j):min(ind0_eq_w(j)+wlen,length(psd));
    i_all0=[i_all0,irange];
end

N7=length(ind1);
ds=eq_dnum(ind1);

%2 x wlen culling for this magnitude range
for i=1:N7
    [~,ind1_store(i,1)]=min(abs(dates_sort-ds(i)));
end
ind1_eq_w=ind1_store(ind1_store>1);

i_all1=[];

for j=1:length(ind1_eq_w)
    irange=ind1_eq_w(j):min(ind1_eq_w(j)+2*wlen,length(psd));
    i_all1=[i_all1,irange];
end

%4 x wlen culling for this magnitude range  
N7_5=length(ind2);
    ds=eq_dnum(ind2);
for i=1:N7_5
    [~,ind2_store(i,1)]=min(abs(dates_sort-ds(i)));
end
ind2_eq_w=ind2_store(ind2_store>1);

i_all2=[];

for j=1:length(ind2_eq_w)
    irange=ind2_eq_w(j):min(ind2_eq_w(j)+4*wlen,length(psd));
    i_all2=[i_all2,irange];
end

%8 x wlen culling for this magnitude range
N8=length(ind3);
    ds=eq_dnum(ind3);
for i=1:N8
    [~,ind3_store(i,1)]=min(abs(dates_sort-ds(i)));
end
ind3_eq_w=ind3_store(ind3_store>1);

i_all3=[];

for j=1:length(ind3_eq_w)
    irange=ind3_eq_w(j):min(ind3_eq_w(j)+8*wlen,length(psd));
    i_all3=[i_all3,irange];
end

%16 x wlen culling for this magnitude range
N9=length(ind4);
    ds=eq_dnum(ind4);
for i=1:N9
    [~,ind4_store(i,1)]=min(abs(dates_sort-ds(i)));
end
ind4_eq_w=ind4_store(ind4_store>1);

i_all4=[];

for j=1:length(ind4_eq_w)
    irange=ind4_eq_w(j):min(ind4_eq_w(j)+16*wlen,length(psd));
    i_all4=[i_all4,irange];
end

%Peterson (1991) noise model culling calculation; performed acceleration integration in same
%period/frequency bands as data

%mid_periods (s)

zlnm = [0.0168    0.0331    0.0679    0.0644    0.0945    0.2425    0.9750 6.3476];
zhnm = [0.0010    0.0033    0.0086    0.0149    0.0240    0.0475    0.2275 1.0884]*1.0e+03;


if vel2
%convert Peterson integrals to velocity if desired
for i=1:8
    zlnm(:,i)=(zlnm(:,i)*(psd_periods(i)/(2*pi)));
    zhnm(:,i)=(zhnm(:,i)*(psd_periods(i)/(2*pi)));
end
end

zl_thresh=sum(zlnm(usebands));
zh_thresh=sum(zhnm(usebands));
n_facl=5;
n_fach=50;

i_all=[i_all0,i_all1,i_all2,i_all3,i_all4];

%remove eq windows flagged above
psd_1=psd;
psd_1(i_all)=NaN;

%Peterson Culling
psd_2=psd_1;
%low power
ind=find(psd_2 < zl_thresh*n_facl);
psd_2(ind)=NaN;

%high power
ind=find(psd_2 > zh_thresh*n_fach);
psd_2(ind)=NaN;

if showfigs
%figure showing all data and earthquakes
figure(1)
semilogy(dates_sort,psd,'k.','markersize',6)
title([st_name,' (',NET,'; ',GLOC,') ',bstr,' s'],'fontsize',40)
hold on
H0=semilogy(dates_sort(i_all0),psd(i_all0),'o','markersize',6,'markeredgecolor',[0.4940 0.1840 0.5560]);
H1=semilogy(dates_sort(i_all1),psd(i_all1),'bo','markersize',7);
H2=semilogy(dates_sort(i_all2),psd(i_all2),'co','markersize',8);
H3=semilogy(dates_sort(i_all3),psd(i_all3),'go','markersize',9);
H4=semilogy(dates_sort(i_all4),psd(i_all4),'ro','markersize',10);
hold off
hold on

%red dot-dash lines on the plot are the peterson limits for culling
plot([min(dates_sort),max(dates_sort)],zl_thresh*n_facl*[1;1],'r-.','linewidth',2)
plot([min(dates_sort),max(dates_sort)],zh_thresh*n_fach*[1;1],'r-.','linewidth',2)

hold off
lh=legend([H0 H1 H2 H3 H4],...
    {'6.0 \leq M_{ww} < 7.0','7.0 \leq M_{ww} < 7.5','7.5 \leq M_{ww} < 8.0','8.0 \leq M_{ww} < 9.0','M_{ww} \geq 9.0'},'fontsize',30,'location','southeast');
set(lh,'color',[0.8 0.8 0.8])
datetick('x')
xlim([min(dates_sort),max(dates_sort)])
ylim([0.001,max(psd_raw)])
ylabel(unitA)
grid on
set(gca,'FontSize',40,'LineWidth',1.0);

end

if showfigs
%Eq and Peterson Culled Hourly Data
figure(2)
semilogy(dates_sort,psd_2,'k.')
datetick('x')
hold off
ylabel(unitA)
title([st_name,' Eq. and Pet. Culled Data (unsmoothed) ',bstr,' s'])
bookfonts
end

if vel2
    psd=psd.^2;
    psd_1=psd_1.^2;
    psd_2=psd_2.^2;
end

%Smoothing Interval in daysfrom raw semi-hourly PSDs for subsequent fitting
%deltad=36.5242/4;

%2 months smoothing for primary on Nat. Commun. paper (~61 days)
deltad=round(T0/6);


%ddy will be the smoothed data timebase
ddy=(min(dates_sort)+deltad/2:1:max(dates_sort)-deltad/2)';

%there are 48 1-hour PSDs/day in the raw data
psd_cull_dyx=movmedian(psd,deltad*48,'omitmissing','Endpoints','shrink');
psd_cull_dy(isnan(psd))=NaN;
psd_cull_dy=interp1(dates_sort,psd_cull_dyx,ddy);

psd_cull_dy_1x=movmedian(psd_1,deltad*48,'omitmissing','Endpoints','shrink');
psd_cull_dy_1(isnan(psd_1))=NaN;
psd_cull_dy_1=interp1(dates_sort,psd_cull_dy_1x,ddy);

psd_cull_dy_2x=movmedian(psd_2,deltad*48,'omitmissing','Endpoints','shrink');
psd_cull_dy_2(isnan(psd_2))=NaN;
psd_cull_dy_2=interp1(dates_sort,psd_cull_dy_2x,ddy);

%conservatively censor out any NaN-related period edge effects by noting where the raw
%data are absent
psd_2_i=psd_2;
psd_2_i(isnan(psd_2))=0;
psdinterp=interp1(dates_sort,psd_2_i,ddy,'nearest');
ind_NaN=find(psdinterp==0);
psd_cull_dy(ind_NaN)=NaN;
psd_cull_dy_1(ind_NaN)=NaN;
psd_cull_dy_2(ind_NaN)=NaN;

%index of all non-NaN data for fitting
ind=find(~isnan(psd_cull_dy));

%robust fitting for variously culled psd data integrals assigned above
[B,B_stat]=robustfit(ddy(ind),psd_cull_dy(ind));

Bstds=[sqrt(B_stat.covb(1,1)),sqrt(B_stat.covb(2,2))];

ind=find(~isnan(psd_cull_dy_1));
[B1,B1_stat]=robustfit(ddy(ind),psd_cull_dy_1(ind));

B1stds=[sqrt(B_stat.covb(1,1)),sqrt(B1_stat.covb(2,2))];

ind=find(~isnan(psd_cull_dy_2));

[B2,B2_stat]=robustfit(ddy(ind),psd_cull_dy_2(ind));

B2stds=[sqrt(B_stat.covb(1,1)),sqrt(B2_stat.covb(2,2))];

nanmeanB=nanmean(psd_cull_dy);
nanmedianB=nanmedian(psd_cull_dy);
pctperyearB=B(2)*T0/nanmedianB*100;
pctperyearB_std=Bstds(2)*T0/nanmedianB*100;

nanmeanB1=nanmean(psd_cull_dy_1);
nanmedianB1=nanmedian(psd_cull_dy_1);
pctperyearB1=B1(2)*T0/nanmedianB1*100;
pctperyearB1_std=B1stds(2)*T0/nanmedianB1*100;

nanmeanB2=nanmean(psd_cull_dy_2);
nanmedianB2=nanmedian(psd_cull_dy_2);

pctperyearB2=B2(2)*T0/nanmedianB2*100;
pctperyearB2_std=B2stds(2)*T0/nanmedianB2*100;

if showfigs
%figure showing successive culling of data and trend fitting (Supplemental
%Fig. 11)
figure(3)
subplot(4,1,1)
plot(ddy,psd_cull_dy,'b.')
grid on
hold on
plot(ddy,B(1)+ddy*B(2),'r','linewidth',3)
hold off
datetick('x')
ylabel(unitA)
set(gca,'FontSize',20,'LineWidth',1.0);

title([num2str(deltad),'-day median; Unculled ',st_name,' ',bstr,' ',num2str(B(2)*T0,2),'\pm'...
    ,num2str(T0*Bstds(2),2),' ',unitR,' (',num2str(pctperyearB,2),'\pm',num2str(pctperyearB_std,2),' %y^{-1})'])
bookfonts

subplot(4,1,2)
plot(ddy,psd_cull_dy_1,'b.')
hold on
plot(ddy,B1(1)+ddy*B1(2),'r','linewidth',3)
grid on
hold off
datetick('x')
ylabel(unitA)

title(['EQ Culled ',st_name,' ',bstr,'s; ',num2str(B1(2)*T0,2),'\pm'...
    ,num2str(T0*B1stds(2),2),' ',unitR,' (',num2str(pctperyearB1,2),'\pm',num2str(pctperyearB1_std,2),' %y^{-1})'])
set(gca,'FontSize',20,'LineWidth',1.0);

subplot(4,1,3)
plot(ddy,psd_cull_dy_2,'b.')
hold on
plot(ddy,B2(1)+ddy*B2(2),'r','linewidth',3)
grid on
hold off
datetick('x')
ylabel(unitA)

title(['EQ and Pet. Culled ',st_name,' ',bstr,'s; ',num2str(B2(2)*T0,2),'\pm'...
    ,num2str(T0*B2stds(2),2),' ',unitR,' (',num2str(pctperyearB2,2),'\pm',num2str(pctperyearB2_std,2),' %y^{-1})'])
set(gca,'FontSize',20,'LineWidth',1.0);

end

%remove annual harmonics (with gaps accomodated, trend removed)
d_base=ddy-datenum(2015,1,1);


%Four annual harmonics are sufficient to remove the stationary component of the annual
%periodicity from the fully culled date series
nharm=4;

%Data (including harmonics) with trend removed
X=psd_cull_dy_2-(B2(1)+ddy*B2(2));

%sine basis functions are aligned with January 1 in phase here via d_base
for i=1:nharm
    sfun(:,i)=sin(2*pi*i*d_base/T0);
    cfun(:,i)=cos(2*pi*i*d_base/T0);
end

not_NaN_ind=find(~isnan(psd_cull_dy_2));

hfun=zeros(size(d_base));

for i=1:nharm
%standard fourier coefficientcalculation
    Norm=length(ind)/2;
    a(i)=dot(sfun(ind,i),X(ind))/Norm;
    b(i)=dot(cfun(ind,i),X(ind))/Norm;


    hfun=hfun+a(i)*sfun(:,i)+b(i)*cfun(:,i);
end
    %insert data gaps into the annual harmonic functional fit
    hfun(~not_NaN_ind)=NaN;

%show harmonic function and data

psd2_deharm=psd_cull_dy_2-hfun;
[cDataSort,csort_ind]=sort(psd2_deharm);
N=length(not_NaN_ind);
indc95=csort_ind(floor(0.95*N):end);

%recalculate slope of the fully culled data with harmonics removed
[B2h,B2h_stat]=robustfit(ddy,psd_cull_dy_2-hfun);
B2hstds=[sqrt(B2h_stat.covb(1,1)),sqrt(B2h_stat.covb(2,2))];

nanmeanB2h=nanmean(psd_cull_dy_2-hfun);
nanmedianB2h=nanmedian(psd_cull_dy_2-hfun);

pctperyearB2h=B2h(2)*T0/nanmedianB2h*100;
pctperyearB2h_std=B2hstds(2)*T0/nanmedianB2h*100;

%recalculate parameters for only 21st century data
ind2000=find(ddy>=datenum(2000,1,1));
if numel(ind2000)>2000
[B2h2000,B2h2000_stat]=robustfit(ddy(ind2000),psd_cull_dy_2(ind2000)-hfun(ind2000));
B2h2000stds=[sqrt(B2h2000_stat.covb(1,1)),sqrt(B2h2000_stat.covb(2,2))];

nanmeanB2h2000=nanmean(psd_cull_dy_2(ind2000)-hfun(ind2000));
nanmedianB2h2000=nanmedian(psd_cull_dy_2(ind2000)-hfun(ind2000));
else
B2h2000=[0 0];
B2h2000stds=[0 0];
nanmedianB2h2000=1;
end

pctperyearB2h2000=B2h2000(2)*T0/nanmedianB2h2000*100;
pctperyearB2h2000_std=B2h2000stds(2)*T0/nanmedianB2h2000*100;

if showfigs
figure(3)
    subplot(4,1,4)
plot(ddy,psd_cull_dy_2-hfun,'b.')
hold on
plot(ddy,B2h(1)+ddy*B2h(2),'r','linewidth',3)
plot(ddy(ind2000),B2h2000(1)+ddy(ind2000)*B2h2000(2),'k--','linewidth',3)
grid on
hold off
axis tight
datetick('x')
title(['EQ and Pet. Culled; H(t) Removed ',st_name,' ',bstr,' ',num2str(B2h(2)*T0,2),'\pm'...
    ,num2str(T0*B2hstds(2),2),' ',unitR,' (',num2str(pctperyearB2h,2),'\pm',num2str(pctperyearB2h_std,2),' %y^{-1})'])
ylabel(unitA)
set(gca,'FontSize',20,'LineWidth',1.0);
dyplot=5*mad(psd_cull_dy_2-hfun);
ylim([nanmedianB2h-dyplot,nanmedianB2h+dyplot])

if showfigs
%Shows time series with harmonic function H(t) and trend
%removed
figure(4)
subplot(4,1,1)
plot(ddy,psd_cull_dy_2,'b.')
hold on
plot(ddy,B2(1)+ddy*B2(2),'r','linewidth',3)
grid on
hold off
datetick('x')
ylabel(unitA)
title([st_name,' ',bstr,'s; ',num2str(B2(2)*T0,2),'\pm'...
    ,num2str(T0*B2stds(2),2),' ',unitR,' (',num2str(pctperyearB2,2),'\pm',num2str(pctperyearB2_std,2),' %y^{-1})'])
set(gca,'FontSize',20,'LineWidth',1.0);
bookfonts

subplot(4,1,2)
plot(ddy(ind),hfun(ind),'.');
datetick('x')
title([num2str(nharm),' Harmonics; Annual Phase: ',num2str(atan2d(a(1),b(1)),2),' Deg.'])
ylabel(unitA)
bookfonts

subplot(4,1,3)
plot(ddy,psd_cull_dy_2-hfun,'b.')
hold on
plot(ddy,B2h(1)+ddy*B2h(2),'r','linewidth',3)
%plot(ddy(ind2000),B2h2000(1)+ddy(ind2000)*B2h2000(2),'k--','linewidth',3)
grid on
hold off
axis tight
datetick('x')
title(['Harmonics Removed; ',num2str(B2h(2)*T0,2),'\pm'...
    ,num2str(T0*B2hstds(2),2),' ',unitR,' (',num2str(pctperyearB2h,2),'\pm',num2str(pctperyearB2h_std,2),' %y^{-1})'])
ylabel(unitA)
set(gca,'FontSize',20,'LineWidth',1.0);
dyplot=5*mad(psd_cull_dy_2-hfun);
ylim([nanmedianB2h-dyplot,nanmedianB2h+dyplot])
bookfonts

subplot(4,1,4)
plot(ddy,psd_cull_dy_2-hfun-(B2h(1)+ddy*B2h(2)),'b.')
hold on
grid on
hold off
axis tight
datetick('x')
title(['Harmonics and Trend Removed; ',num2str(B2h(2)*T0,2),'\pm'...
    ,num2str(T0*B2hstds(2),2),' ',unitR,' (',num2str(pctperyearB2h,2),'\pm',num2str(pctperyearB2h_std,2),' %y^{-1})'])
ylabel(unitA)
set(gca,'FontSize',20,'LineWidth',1.0);
dyplot=5*mad(psd_cull_dy_2-hfun);
ylim([-dyplot,dyplot])
bookfonts
set (gcf,'position',[924 219 824 1143]);
end

end

%end master station loop below

%number of years for output text
Nyears=(max(ddy)-min(ddy))/T0;

%%%%%

%completeness fraction for output text
Cmp_frac=numel(psd_cull_dy_2(~isnan(psd_cull_dy_2)))/numel(psd_cull_dy_2);

if B2h2000(1)==0
    Cmp_frac2000=0;
else
Cmp_frac2000=numel(psd_cull_dy_2(~isnan(psd_cull_dy_2(ind2000))))/numel(psd_cull_dy_2(ind2000));
end
disp([num2str(ista),' ',st_name,' (',NET,'; ',GLOC,') ',num2str(Nyears,3),' ',num2str(Cmp_frac),' (',num2str(Cmp_frac2000),')'])
disp(['EQ and Pet. Culled; Ann. Harmonics Removed ',st_name,' ',bstr,' ',num2str(B2h(2)*T0,2),' +- '...
    ,num2str(T0*B2hstds(2),2),' ',unitR,' (',num2str(pctperyearB2h,2),'+-',num2str(pctperyearB2h_std,2),' %y^{-1})'])
disp(['EQ and Pet. Culled; Ann. Harmonics Removed (post 2000) ',st_name,' ',bstr,' ',num2str(B2h2000(2)*T0,2),' +- '...
    ,num2str(T0*B2h2000stds(2),2),' ',unitR,' (',num2str(pctperyearB2h2000,2),'+-',num2str(pctperyearB2h2000_std,2),' %y^{-1})'])

%Save all results as below in the Results directory for subsequent analyses
if acc
    eval(['save Results/results2.0_',bstrf,st_name,'_acc.mat st_name latlon T0 B B1 B2 B2h* Bstds B1stds B2stds B2hstds* Cmp_* nanme* pctperyear* psd_cull_dy* hfun ddy deltad dates_sort psd_2 n_fac* a b'])
end

if vel2
    eval(['save Results/results2.0_',bstrf,st_name,'_vel2.mat st_name latlon T0 B B1 B2 B2h* Bstds B1stds B2stds B2hstds* Cmp_* nanme* pctperyear* psd_cull_dy* hfun ddy deltad dates_sort psd_2 n_fac* a b'])
end

end




