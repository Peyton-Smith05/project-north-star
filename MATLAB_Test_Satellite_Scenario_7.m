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


% Visualize the Satellites and their Orbits
v = satelliteScenarioViewer(sc);

% Visualize a Dynamic Animation of the Satellite Movement
play(sc)


% Obtain the Position and Velocity History of the Satellites
[positionSGP4, velocitySGP4, time] = states(satSGP4);
positionSGP4_magnitude = vecnorm(positionSGP4,2,1);
velocitySGP4_magnitude = vecnorm(velocitySGP4,2,1);


% % Plot Plots of Position and Velocity History of Satellites
% % Individual Plots
% for index = (1:number_of_satellites)
%     plot(time, positionSGP4_magnitude(:,:,index))
%     xlabel("Time")
%     ylabel("Position (m)")
%     legend(string(index) + " SGP4")
%     saveas(gcf, "plot_position_" + string(index) + ".png")
%     
%     plot(time, velocitySGP4_magnitude(:,:,index))
%     xlabel("Time")
%     ylabel("Velocity deviation (m/s)")
%     legend(string(index) + " SGP4")
%     saveas(gcf, "plot_velocity_" + string(index) + ".png")
% end
% 
% % Group Plots
% for index = (1:number_of_satellites)
%     plot(time, positionSGP4_magnitude(:,:,index))
%     xlabel("Time")
%     ylabel("Position (m)")
%     hold on
% end
% % legend(string(1:number_of_satellites) + " SGP4")
% saveas(gcf, "plot_position.png")
% close all
% 
% for index = (1:number_of_satellites)
%     plot(time, velocitySGP4_magnitude(:,:,index))
%     xlabel("Time")
%     ylabel("Velocity deviation (m/s)")
%     hold on
% end
% % legend(string(1:number_of_satellites) + " SGP4")
% saveas(gcf, "plot_velocity.png")
% close all

