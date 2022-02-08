
function [data] = windsurfer(targetLat,targetLon,targetYear,targetMonth,...
    targetDay,output,resolution,varargin)
%
% [data] = windsurfer(targetLat,targetLon,targetYear,targetMonth,...
%     targetDay,targetHour,output,resolution)
% downloads and processes automatically the WINDSURFER/NORA3 wave field hindcast with a  time
% resolution of 1 h and a spatial resolution of ca. 3 km.
% The data set covers the period from 1993 to 2020. The wave hindcast are
% obtained with WAM wave model version cycle 4.7.0.
% The Norwegian Meteorological Institute provided the data freely.
% These are available on their Thredds server [1].
% The documentation is contained in the NetCDF files that are read
% through the OPeNDAP  framework.
% [1] https://thredds.met.no/thredds
%
% Input:
%    * targetLat: [1x1] or [1xNlat] double: target latitude where the data
%         should be extracted
%    * targetLon: [1x1] or [1xNlon] double: target longitude where the data
%         should be extracted
%    * targetYear: [1x1] or [1xNyear] double: year(s) from which the data are
%         taken
%    * targetMonth: [1x1] or [1xNmonth] double: month(s) from which the data are
%         taken
%    * targetDay: [1x1] or [1xNday] double: day(s) from which the data are
%         taken
%    * targetHour: [1x1] or [1xNhour] double: hour(s) from which the data are
%         taken
%    * output: cell of string: Name of the variables to read and extract
%    from the netcdf file
%    resolution: [1x1]:double: resolution, in degree for the gridded data
%
% Outputs:
%   * data: structure with the following fields
%     - 'time': [1x1] datetime
%     - 'lon': [Nlon x Nlat] double: longitude (deg)
%     - 'lat': [Nlon x Nlat] double: latitude (deg)
%     - 'ff': [Nlon x Nlat] double: wind speed at 10 m asl
%     - 'dd': [Nlon x Nlat] double: wind direction at 10 m asl
%     - 'fv': [Nlon x Nlat] double: friction velocity (m/s)
%     - 'dc': [Nlon x Nlat] double: drag coefficient ()
%     - 'hs': [Nlon x Nlat] double:  Total significant wave height
%     - 'tp': [Nlon x Nlat] double:  Total peak period
%     - 'hs_sea': [Nlon x Nlat] double: sea surface wind wave significant
%     height (m)
%     - 'tp_sea': [Nlon x Nlat] double: sea surface wind wave peak period
%     from variance spectral density (s)
%     - 'hs_swell': [Nlon x Nlat] double:
%     sea_surface swell wave significant height (m)
%     - 'tp_swell': [Nlon x Nlat] double:
%     sea surface swell wave peak period from variance spectral density (s)
%     - 'thq': [Nlon x Nlat] double: Total mean wave direction (deg)
%     - 'thq_sea': [Nlon x Nlat] double: wind wave mean wave direction (deg)
%     - 'thq_swell:' Swell wave mean wave direction (deg)
%     - 'sic': [Nlon x Nlat] double:  sea ice area fraction
%     - 'sit': [Nlon x Nlat] double: sea ice thickness (m)
%     - 'model_depth': [Nlon x Nlat] double:
%     sea floor depth below sea level (m)
%
% Author: E. Cheynet - UiB, Norway - last modified: 06-12-2021

%% Optional aprameters


p = inputParser();
p.CaseSensitive = false;
p.addOptional('newLon',[]); % optional parameters -- longitude vector
p.addOptional('newLat',[]); % optional parameters -- latitude vector
p.parse(varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%
newLat = p.Results.newLat ;
newLon = p.Results.newLon ;

%% Definition of the grid
if ~isempty(newLon) && ~isempty(newLat)
    lon = newLon;
    lat = newLat;
    Nlat = numel(lat);
else
    [lon,lat] = meshgrid(targetLon(1):resolution:targetLon(end),targetLat(1):resolution:targetLat(end));
    [Nlat,Nlon]=size(lat);
end

if numel(targetDay)==1, targetDay = ['0',targetDay];end
if numel(targetMonth)==1, targetMonth = ['0',targetMonth];end
% if numel(targetHour)==1, targetHour = ['0',targetHour];end
%% Preallocation and initalisation
data = struct('time',[],'lon',[],'lat',[],'ff',[],'dd',[],'fv',[],'dc',[],'hs',[],'tp',[],...
    'hs_sea',[],'tp_sea',[],'hs_swell',[],'tp_swell',[],...
    'thq',[],'thq_sea',[],'thq_swell',[],'sic',[],'sit',[],'model_depth',[]);
urldat= ['https://thredds.met.no/thredds/dodsC/windsurfer/mywavewam3km_files/',...
    targetYear,'/',targetMonth,'/',targetYear,targetMonth,targetDay,'_MyWam3km_hindcast.nc'];
time0 = ncread(urldat,'time')./86400+datenum('1970-01-01 00:00:00');
time = datetime(datestr(double(time0)));
% [~,indTime]=min(abs(hour(time)-str2double(targetHour)));
data.time = time;
%% Get narrower area to speed up the reading procedure
lon00 = ncread(urldat,'longitude');
lat00 = ncread(urldat,'latitude');


if numel(targetLat)==1 &&   numel(targetLon)==1
    offsetLat = 0.3;
    offsetLon = 0.15;
    [~,indStart] = min(sqrt((lat00(:)-(targetLat(1)-offsetLat)).^2 + abs(lon00(:)-(targetLon(1)-offsetLon)).^2));
    [~,indEnd] = min(sqrt((lat00(:)-(targetLat(1)+offsetLat)).^2 + abs(lon00(:)-(targetLon(1)+offsetLon)).^2));
elseif numel(targetLat)==2 &&   numel(targetLat)==2
    
    offsetLat = 0.5*sqrt(diff(targetLat).^2 + diff(targetLon).^2);
    offsetLon = offsetLat;
    [~,indStart] = min(sqrt((lat00(:)-(targetLat(1)-offsetLat)).^2 + abs(lon00(:)-(targetLon(2)-offsetLon)).^2));
    [~,indEnd] = min(sqrt((lat00(:)-(targetLat(2)+offsetLat)).^2 + abs(lon00(:)-(targetLon(1)+offsetLon)).^2));
else
    error('targetLat and targetLatmust have the same dimensions')
end


[row1, col1] = ind2sub(size(lat00), indStart);
[row2, col2] = ind2sub(size(lat00), indEnd);
r1 = min(row1,row2); % row start
cr = abs(row2-row1+1);% row count
c1 = min(col1,col2); % column start
cc = abs(col2-col1+1); % column count
lon00 = ncread(urldat,'longitude',[r1,c1],[cr,cc]);
lat00 = ncread(urldat,'latitude',[r1,c1],[cr,cc]);
dummyLat = double(lat00(:));
dummyLon = double(lon00(:));
ind = find(dummyLat>=min(targetLat(:)-offsetLat) & dummyLat <= max(targetLat(:)+offsetLat) &...
    dummyLon>=min(targetLon(:)-offsetLon) & dummyLon <= max(targetLon(:)+offsetLon));


%% Read the data in a for loop for each selected output
%  Data are resampled spatially as gridded data

    
Nout = numel(output);
Nhour = 24; % number of hours for 1 day

if ~isempty(newLat),
    for ii=1:Nout
        if strcmpi(output{ii},'model_depth'),
            myVar0 = ncread(urldat,output{ii},[r1,c1],[cr,cc]);
            dummyVar = double(reshape(myVar0,[],1));
            F_Var = scatteredInterpolant(dummyLat(ind),dummyLon(ind),dummyVar(ind),'linear','none');
            data.(lower(output{ii}))(jj,:) = F_Var(lat,lon);
        else
            myVar0 = ncread(urldat,output{ii},[r1,c1,1],[cr,cc,Nhour]);
            
            
            for jj=1:Nhour
                dummyVar = double(reshape(myVar0(:,:,jj),[],1));
                if strcmpi(output{ii},'dd') || contains(output{ii},'thq')
                    un = cosd(dummyVar);
                    ue = sind(dummyVar);
                    F_ue = scatteredInterpolant(dummyLat(ind),dummyLon(ind),un(ind),'linear','none');
                    F_un = scatteredInterpolant(dummyLat(ind),dummyLon(ind),ue(ind),'linear','none');
                    newUe = F_ue(lat,lon);
                    newUn = F_un(lat,lon);
                    newDir = atan2(newUn,newUe).*180/pi;
                    newDir(newDir<=0) = newDir(newDir<=0)+ 360;
                    data.(lower(output{ii}))(jj,:) = newDir;
                else
                    F_Var = scatteredInterpolant(dummyLat(ind),dummyLon(ind),dummyVar(ind),'linear','none');
                    data.(lower(output{ii}))(jj,:) = F_Var(lat,lon);
                end
            end
        end
    end
else
    for ii=1:Nout
        if strcmpi(output{ii},'model_depth'),
            myVar0 = ncread(urldat,output{ii},[r1,c1],[cr,cc]);
            dummyVar = double(reshape(myVar0,[],1));
            F_Var = scatteredInterpolant(dummyLat(ind),dummyLon(ind),dummyVar(ind),'linear','none');
            data.(lower(output{ii}))(jj,:,:) = F_Var(lat,lon);
        else
            myVar0 = ncread(urldat,output{ii},[r1,c1,1],[cr,cc,Nhour]);
            
            
            for jj=1:Nhour
                dummyVar = double(reshape(myVar0(:,:,jj),[],1));
                if strcmpi(output{ii},'dd') || contains(output{ii},'thq')
                    un = cosd(dummyVar);
                    ue = sind(dummyVar);
                    F_ue = scatteredInterpolant(dummyLat(ind),dummyLon(ind),un(ind),'linear','none');
                    F_un = scatteredInterpolant(dummyLat(ind),dummyLon(ind),ue(ind),'linear','none');
                    newUe = F_ue(lat,lon);
                    newUn = F_un(lat,lon);
                    newDir = atan2(newUn,newUe).*180/pi;
                    newDir(newDir<=0) = newDir(newDir<=0)+ 360;
                    data.(lower(output{ii}))(jj,:,:) = newDir;
                else
                    F_Var = scatteredInterpolant(dummyLat(ind),dummyLon(ind),dummyVar(ind),'linear','none');
                    data.(lower(output{ii}))(jj,:,:) = F_Var(lat,lon);
                end
            end
        end
    end
end
data.lon = lon;
data.lat = lat;

end
