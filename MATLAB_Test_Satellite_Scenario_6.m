% Create Satellite Scenario
startTime = datetime(2022,5,11,12,35,38);
stopTime = startTime + days(2); % 48 hours length
sampleTime = 60; % seconds (I think)
sc = satelliteScenario(startTime,stopTime,sampleTime)


% Add Satellites to the Satellite Scenario
tleFile = "leoSatelliteConstellation.tle";
lines_per_satellite = 3 % TLE files usually have a new satellite every 3 lines
number_of_satellites = length(readlines(tleFile)) / lines_per_satellite
Name = string(1:number_of_satellites);
SGP4_Name = Name + 'SGP4';


satSGP4 = satellite(sc,tleFile, ...
    "Name",SGP4_Name, ...
    "OrbitPropagator","sgp4") % Red


% % Visualize the Satellites and their Orbits
% v = satelliteScenarioViewer(sc);
% 
% % Visualize a Dynamic Animation of the Satellite Movement
% play(sc)


% Obtain the Position and Velocity History of the Satellites
[positionSGP4, velocitySGP4, time] = states(satSGP4);
positionSGP4_magnitude = vecnorm(positionSGP4,2,1);
velocitySGP4_magnitude = vecnorm(velocitySGP4,2,1);

for index = (1:number_of_satellites)
    plot(time, positionSGP4_magnitude(:,:,index))
    xlabel("Time")
    ylabel("Position (m)")
    legend(string(index) + " SGP4")
    saveas(gcf, "plot_position_" + string(index) + ".png")
    
    plot(time, velocitySGP4_magnitude(:,:,index))
    xlabel("Time")
    ylabel("Velocity deviation (m/s)")
    legend(string(index) + " SGP4")
    saveas(gcf, "plot_velocity_" + string(index) + ".png")
end

% test_index = 40
% plot(time, positionSGP4_magnitude(:,:,test_index))
% xlabel("Time")
% ylabel("Position (m)")
% legend(string(test_index) + " SGP4")
% saveas(gcf, "plot_position_" + string(test_index) + ".png")
% 
% plot(time, velocitySGP4_magnitude(:,:,test_index))
% xlabel("Time")
% ylabel("Velocity deviation (m/s)")
% legend(string(test_index) + " SGP4")
% saveas(gcf, "plot_velocity_" + string(test_index) + ".png")


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