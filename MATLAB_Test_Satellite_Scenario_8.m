% Create Satellite Scenario
startTime = datetime(2022,6,11,12,35,38);
stopTime = startTime + days(3); % hours length
sampleTime = 60; % seconds between refresh
sc = satelliteScenario(startTime,stopTime,sampleTime)


% Add Satellites to the Satellite Scenario
tleFile = "22-144 LEO ELSETs Copy.txt";
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

% % For showing the number of possible conjunctions
% count
% possible_conj_indeces


% Verify Possible Conjunctions
count_2 = 0;
coordinate_diff = [];
time_stamp = 0;
satellite_A_number = 0;
satellite_A_coordinates = [];
satellite_B_number = 0;
satellite_B_coordinates = [];
actual_conjunctions = []; % satellite#A, satellite#B, timestamp
near_coordinates = near / sqrt(3);

for index_6 = (1:length(possible_conj_indeces(:,3))) % Length of each timestamp
    
    time_stamp = possible_conj_indeces(index_6,3); % current timestamp

    satellite_A_number = possible_conj_indeces(index_6,1);
    satellite_A_coordinates = ...
        pos_by_timepages(satellite_A_number,:,time_stamp);
    
    satellite_B_number = possible_conj_indeces(index_6,2);
    satellite_B_coordinates = ...
        pos_by_timepages(satellite_B_number,:,time_stamp);
    
    coordinate_diff = satellite_A_coordinates - satellite_B_coordinates;
    coordinate_diff = abs(coordinate_diff); % turn them all to absolute values to compare

    if (coordinate_diff(1,1) < near_coordinates) ...
        && (coordinate_diff(1,2) < near_coordinates) ...
        && (coordinate_diff(1,3) < near_coordinates)
        
        count_2 = count_2 + 1;
        actual_conjunctions(count_2,:) = possible_conj_indeces(index_6,:);
    end
end

count_2
actual_conjunctions