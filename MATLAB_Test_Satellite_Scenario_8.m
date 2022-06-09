% Create Satellite Scenario
startTime = datetime(2022,5,11,12,35,38);
stopTime = startTime + hours(21/60); % hours length
sampleTime = 60; % seconds between refresh
sc = satelliteScenario(startTime,stopTime,sampleTime)


% Add Satellites to the Satellite Scenario
tleFile = "leoSatelliteConstellation2.tle";
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

trans_time = transpose(time); % time as a column
trans_position = pagetranspose(positionSGP4); % by time column with one satellite per page
trans_pos_mag = pagetranspose(positionSGP4_magnitude); % by time column with one satellite per page

% arrays of position magnitudes by satellite column with one time stamp per
% page
for index_6 = (1:length(time))
    for index_7 = (1:number_of_satellites)
        pos_mag_by_timepages(index_7,:,index_6) = trans_pos_mag(index_6,:,index_7);
        pos_by_timepages(index_7,:,index_6) = trans_position(index_6,:,index_7);
    end
end
pos_mag_by_timepages;
pos_by_timepages;

% sort magnitudes by 
[mag_fill, mag_indeces] = sort(pos_mag_by_timepages);

% Distance Check
near = 1000 % meters
differences = diff(mag_fill);
% mag_fill(:,:,end);
% differences(:,:,end);

% Find possible conjunctions
count = 0;
possible_conj_indeces = []

for index_4 = (1:length(time)) % Each page is time
    for index_5 = (1:length(differences(:,:,index_4))) % Each row is a satallite
        if differences(index_5,:,index_4) <= near % check difference in magnitude
            count = count + 1;
            possible_conj_indeces(count,1) = mag_indeces(index_5,:,index_4);
            possible_conj_indeces(count,2) = mag_indeces(index_5+1,:,index_4);
            possible_conj_indeces(count,3) = index_4;
        end
    end
end 

count
possible_conj_indeces


% Verify Possible Conjunctions