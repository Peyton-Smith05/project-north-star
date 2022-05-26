% Create Satellite Scenario
startTime = datetime(2022,5,11,12,35,38);
stopTime = startTime + days(2); % 48 hours length
sampleTime = 60;
sc = satelliteScenario(startTime,stopTime,sampleTime)


% Add Satellites to the Satellite Scenario
tleFile = "leoSatelliteConstellation.tle";
m = 3
n = length(readlines(tleFile))/m
Name = string(1:n);
tbk_Name = Name + 'tbk';
SGP4_Name = Name + 'SGP4';
SDP4_Name = Name + 'SDP4';

satTwoBodyKeplerian = satellite(sc,tleFile, ...
    "Name",tbk_Name, ...
    "OrbitPropagator","two-body-keplerian") % Red

satSGP4 = satellite(sc,tleFile, ...
    "Name",SGP4_Name, ...
    "OrbitPropagator","sgp4") % Green

satSDP4 = satellite(sc,tleFile, ...
    "Name",SDP4_Name, ...
    "OrbitPropagator","sdp4") % Pink


% Visualize the Satellites and their Orbits
v = satelliteScenarioViewer(sc);

for i = (1:n) 
    satSGP4(i).MarkerColor = [0 1 0];
    satSGP4(i).Orbit.LineColor = [0 1 0];
    satSGP4(i).LabelFontColor = [0 1 0];
    satSDP4(i).MarkerColor = [1 0 1];
    satSDP4(i).Orbit.LineColor = [1 0 1];
    satSDP4(i).LabelFontColor = [1 0 1];
end

% satSGP4.MarkerColor = [0 1 0];
% satSGP4.Orbit.LineColor = [0 1 0];
% satSGP4.LabelFontColor = [0 1 0];
% satSDP4.MarkerColor = [1 0 1];
% satSDP4.Orbit.LineColor = [1 0 1];
% satSDP4.LabelFontColor = [1 0 1];

% camtarget(v,satSGP4);

% Visualize a Dynamic Animation of the Satellite Movement
play(sc)

% camtarget(v,satSGP4);


% Obtain the Position and Velocity History of the Satellites
[positionTwoBodyKeplerian,velocityTwoBodyKeplerian,time] = states(satTwoBodyKeplerian);
[positionSGP4,velocitySGP4] = states(satSGP4);
[positionSDP4,velocitySDP4] = states(satSDP4);

% 
% % Plot the Magnitude of Relative Position with Respect to Two-Body-Kelerian
% % Prediction
% sgp4RelativePosition = vecnorm(positionSGP4 - positionTwoBodyKeplerian,2,1);
% sdp4RelativePosition = vecnorm(positionSDP4 - positionTwoBodyKeplerian,2,1);
% 
% sgp4RelativePositionKm = sgp4RelativePosition/1000;
% sdp4RelativePositionKm = sdp4RelativePosition/1000;
% plot(time,sgp4RelativePositionKm,time,sdp4RelativePositionKm)
% xlabel("Time")
% ylabel("Relative position (km)")
% legend("SGP4","SDP4")
% saveas(gcf, 'Position_Plot.png')
% 
% % Plot Magnitude of Relative Velocity with Respect to Two-Body-Keplerian
% % Prediction
% sgp4RelativeVelocity = vecnorm(velocitySGP4 - velocityTwoBodyKeplerian,2,1);
% sdp4RelativeVelocity = vecnorm(velocitySDP4 - velocityTwoBodyKeplerian,2,1);
% 
% plot(time,sgp4RelativeVelocity,time,sdp4RelativeVelocity)
% xlabel("Time")
% ylabel("Velocity deviation (m/s)")
% legend("SGP4","SDP4")
% saveas(gcf, 'Velocity_Plot.png')