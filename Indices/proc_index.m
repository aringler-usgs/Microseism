%Aster et al. Nature Communications 2023 Primary Microseism Study
%Processes climate index files for use in clustering
%Time series were downloaded from https://psl.noaa.gov/data/climateindices/list
%
clear
i_names={'BEST','NAO','PDO','PNA','SOI','WP','AMO'};
for j=1:length(i_names)
    i_name=char(i_names(j))
X=importdata([i_name,'.txt']);
[a,b]=size(X);
k=1;
for i=1:a
    for l=2:b
        eval([i_name,'_year(k,1)=X(i)+(l-1.5)/12;']);
        eval([i_name,'_date(k,1)=datenum(',i_name,'_year(k,1),1,1);'])
        eval([i_name,'_data(k,1)=X(i,l);']);
        k=k+1;
    end
end
%interpolate for intercomparison with seismic data, and save as .mat files
dmin=datenum(1985,1,1);
dmax=datenum(2022,8,1);
dinterp=dmin:dmax;
%Apply 61-day smoothing here to climate index data interpolated to daily samples between dmin and dmax
eval([i_name,'_data_interp=movmean(interp1(',i_name,'_date,',i_name,'_data,dinterp,''spline''),61);'])
eval(['save ',i_name,'.mat ',i_name,'_year ',i_name,'_data ',i_name,'_data_interp dinterp'])
end
figure(2000)
plot(dinterp,NAO_data_interp)
datetick('x')
