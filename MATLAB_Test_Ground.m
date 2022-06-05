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
cam = conicalSensor(sat,"Name",names,"MaxViewAngle",90)

% Define the Geographical Site to be Photographed in the Satellite Scenario
name = "Geographical Site";
minElevationAngle = 30; % degrees
geoSite = groundStation(sc, ...
    "Name",name, ...
    "MinElevationAngle",minElevationAngle)

% Add Access Analysis Between the Cameras and the Geographical Site
ac = access(cam,geoSite);
% Properties of access analysis objects
ac(1)

% Visualize the Scenario
v = satelliteScenarioViewer(sc,"ShowDetails",false);
sat(4).ShowLabel = true;
geoSite.ShowLabel = true;
show(sat(4));

% Visualize the Field Of View of the Camera
fov = fieldOfView(cam([cam.Name] == "Satellite 4 Camera"))

% Customize the Visualizations
ac.LineColor = 'green';

% Determine the Times when the Cameras can Photograph the Geographical Site
accessIntervals(ac)
play(sc);

% Calculate System-Wide Access Percentage
