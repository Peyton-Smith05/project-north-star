% Create Satellite Scenario
startTime = datetime(2020,5,12,13,0,0);
stopTime = startTime + hours(6);
sampleTime = 30; % seconds
sc = satelliteScenario(startTime,stopTime,sampleTime)

% Add Satellites to the Satellite Scenario
tleFile = "leoSatelliteConstellation.tle";
sat = satellite(sc,tleFile)

% Add Cameras to the Satellites
names = sat.Name + " Camera";
max_view_angle = 90; % degrees
cam = conicalSensor(sat,"Name",names,"MaxViewAngle",max_view_angle)

% Define the Geographical Site to be Photographed in the Satellite Scenario
latitudes = [42.3001 30];
longitudes = [-71.3504 -100];
num_ground = length(latitudes)
name = string((1:num_ground)) + "Geographical Site";
minElevationAngle = 30; % degrees
geoSite = groundStation(sc,latitudes,longitudes,"Name",name, ...
    "MinElevationAngle",minElevationAngle)

% Add Access Analysis Between the Cameras and the Geographical Site
ac = access(cam,geoSite(2));
% Properties of access analysis objects
ac(1)

% Visualize the Scenario
sat_number = 5;
v = satelliteScenarioViewer(sc,"ShowDetails",false);
sat.ShowLabel = true;
geoSite.ShowLabel = true;
show(sat);

% Visualize the Field Of View of the Camera
fov = fieldOfView(cam())

% Customize the Visualizations
ac.LineColor = 'green';

% Determine the Times when the Cameras can Photograph the Geographical Site
accessIntervals(ac)
play(sc);



% 
% % Calculate System-Wide Access Percentage
% for idx = 1:numel(ac)
%     [s,time] = accessStatus(ac(idx));
%     
%     if idx == 1
%         % Initialize system-wide access status vector in the first iteration
%         systemWideAccessStatus = s;
%     else
%         % Update system-wide access status vector by performing a logical OR
%         % with access status for the current camera-site access
%         % analysis
%         systemWideAccessStatus = or(systemWideAccessStatus,s);
%     end
% end
% 
% plot(time,systemWideAccessStatus,"LineWidth",2);
% grid on;
% xlabel("Time");
% ylabel("System-Wide Access Status");
% 
% n = nnz(systemWideAccessStatus)
% 
% systemWideAccessDuration = n*sc.SampleTime % seconds
% 
% scenarioDuration = seconds(sc.StopTime - sc.StartTime)
% 
% systemWideAccessPercentage = (systemWideAccessDuration/scenarioDuration)*100
% 
% % Improve the System-Wide Access Percentage by Making the Cameras Track the
% % Geographical Site
% pointAt(sat,geoSite);
% 
% % Calculate system-wide access status
% for idx = 1:numel(ac)
%     [s,time] = accessStatus(ac(idx));
%     
%     if idx == 1
%         % Initialize system-wide access status vector in the first iteration
%         systemWideAccessStatus = s;
%     else
%         % Update system-wide access status vector by performing a logical OR
%         % with access status for the current camera-site combination
%         systemWideAccessStatus = or(systemWideAccessStatus,s);
%     end
% end
% 
% % Calculate system-wide access percentage
% n = nnz(systemWideAccessStatus);
% systemWideAccessDuration = n*sc.SampleTime;
% systemWideAccessPercentageWithTracking = (systemWideAccessDuration/scenarioDuration)*100
% 
% play(sc)
