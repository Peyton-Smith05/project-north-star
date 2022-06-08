% Define Mission Parameters and Constellation Initial Conditions
mission.StartDate = datetime(2022, 5, 30, 22, 23, 24);
mission.Duration  = hours(24);

mission.Satellites.SemiMajorAxis  = 29599.8e3 * ones(24,1); % meters
mission.Satellites.Eccentricity   = 0.0005    * ones(24,1);
mission.Satellites.Inclination    = 56        * ones(24,1); % deg
mission.Satellites.ArgOfPeriapsis = 350       * ones(24,1); % deg

mission.Satellites.RAAN               = sort(repmat([0 120 240], 1,8))'; % deg
mission.Satellites.TrueAnomaly        = repmat(0:45:315, 1,3)'; % deg

mission.Satellites.TrueAnomaly(9:16)  = mission.Satellites.TrueAnomaly(9:16)  + 15;
mission.Satellites.TrueAnomaly(17:24) = mission.Satellites.TrueAnomaly(17:24) + 30;

ConstellationDefinition = table(mission.Satellites.SemiMajorAxis, ...
    mission.Satellites.Eccentricity, ...
    mission.Satellites.Inclination, ...
    mission.Satellites.RAAN, ...
    mission.Satellites.ArgOfPeriapsis, ...
    mission.Satellites.TrueAnomaly, ...
    'VariableNames', ["a (m)", "e", "i (deg)", "Ω (deg)", "ω (deg)", "ν (deg)"])


% Open and Configure the Orbit Propagation Model
mission.mdl = "OrbitPropagatorBlockExampleModel";
open_system(mission.mdl);

mission.Satellites.blk = mission.mdl + "/Orbit Propagator";

set_param(mission.Satellites.blk, ...
    "startDate",      num2str(juliandate(mission.StartDate)), ...
    "stateFormatNum", "Orbital elements", ...
    "orbitType",      "Keplerian", ...
    "semiMajorAxis",  "mission.Satellites.SemiMajorAxis", ...
    "eccentricity",   "mission.Satellites.Eccentricity", ...
    "inclination",    "mission.Satellites.Inclination", ...
    "raan",           "mission.Satellites.RAAN", ...
    "argPeriapsis",   "mission.Satellites.ArgOfPeriapsis", ...
    "trueAnomaly",    "mission.Satellites.TrueAnomaly");

set_param(mission.Satellites.blk, ...
    "centralBody",  "Earth", ...
    "outportFrame", "Fixed-frame");

set_param(mission.Satellites.blk, ...
    "propagator",   "Numerical (high precision)", ...
    "gravityModel", "Oblate ellipsoid (J2)", ...
    "useEOPs",      "off");

set_param(mission.mdl, ...
    "SolverType", "Variable-step", ...
    "SolverName", "VariableStepAuto", ...
    "RelTol",     "1e-6", ...
    "AbsTol",     "1e-7", ...
    "StopTime",   string(seconds(mission.Duration)));

set_param(mission.mdl, ...
    "SaveOutput", "on", ...
    "OutputSaveName", "yout", ...
    "SaveFormat", "Dataset");


% Run the Model and Collect Satellite Ephemerides
mission.SimOutput = sim(mission.mdl);

mission.Satellites.TimeseriesPosECEF = mission.SimOutput.yout{1}.Values;
mission.Satellites.TimeseriesVelECEF = mission.SimOutput.yout{2}.Values;

mission.Satellites.TimeseriesPosECEF.TimeInfo.StartDate = mission.StartDate;
mission.Satellites.TimeseriesVelECEF.TimeInfo.StartDate = mission.StartDate;

mission.Satellites.TimeseriesPosECEF


% Load Satellite Ephemerides into satelliteScenario Object
scenario = satelliteScenario(mission.StartDate, mission.StartDate + hours(24), 60);

sat = satellite(scenario, mission.Satellites.TimeseriesPosECEF, mission.Satellites.TimeseriesVelECEF, ...
    "CoordinateFrame", "ecef", "Name", "GALILEO " + (1:24))

disp(scenario)


% Set Graphical Properties on Satellites
set(sat, "ShowLabel", false);
hide([sat(:).GroundTrack]);

set(sat(1:8), "MarkerColor", "red");
set(sat(9:16), "MarkerColor", "blue");
set(sat(17:24), "MarkerColor", "green");
orbit = [sat(:).Orbit];
set(orbit(1:8), "LineColor", "red");
set(orbit(9:16), "LineColor", "blue");
set(orbit(17:24), "LineColor", "green");


% Add Ground Stations to Scenario
gsUS = groundStation(scenario, 42.30048, -71.34908, ...
    "MinElevationAngle", 10, "Name", "Natick");
gsDE = groundStation(scenario, 48.23206, 11.68445, ...
    "MinElevationAngle", 10, "Name", "Munchen");
gsIN = groundStation(scenario, 12.94448, 77.69256, ...
    "MinElevationAngle", 10, "Name", "Bangalore");

figure
geoscatter([gsUS.Latitude gsDE.Latitude gsIN.Latitude], ...
    [gsUS.Longitude gsDE.Longitude gsIN.Longitude], "red", "filled")
geolimits([-75 75], [-180 180])
title("Ground Stations")


% Compute Ground Station to Satellite Access (Line-of-Sight Visibility)
for idx = 1:numel(sat)
    access(gsUS, sat(idx));
    access(gsDE, sat(idx));
    access(gsIN, sat(idx));
end
accessUS = [gsUS(:).Accesses];
accessDE = [gsDE(:).Accesses];
accessIN = [gsIN(:).Accesses];

set(accessUS(1:8), "LineColor", "red");
set(accessUS(9:16), "LineColor", "blue");
set(accessUS(17:24), "LineColor", "green");

set(accessDE(1:8), "LineColor", "red");
set(accessDE(9:16), "LineColor", "blue");
set(accessDE(17:24), "LineColor", "green");

set(accessIN(1:8), "LineColor", "red");
set(accessIN(9:16), "LineColor", "blue");
set(accessIN(17:24), "LineColor", "green");

intervalsUS = accessIntervals(accessUS);
intervalsUS = sortrows(intervalsUS, "StartTime", "ascend")

intervalsDE = accessIntervals(accessDE);
intervalsDE = sortrows(intervalsDE, "StartTime", "ascend")

intervalsIN = accessIntervals(accessIN);
intervalsIN = sortrows(intervalsIN, "StartTime", "ascend")


% View the Satellite Scenario
viewer3D = satelliteScenarioViewer(scenario);


% Compare Access Between Ground Station
% Initialize array with size equal to number of timesteps in scenario
timeSteps = mission.StartDate:seconds(60):mission.StartDate+days(1);
statusUS = zeros(1, numel(timeSteps));
statusDE = statusUS;
statusIN = statusUS;

% Sum cumulative access at each timestep
for idx = 1:24
    statusUS = statusUS + accessStatus(accessUS(idx));
    statusDE = statusDE + accessStatus(accessDE(idx));
    statusIN = statusIN + accessStatus(accessIN(idx));
end
clear idx;

subplot(3,1,1);
plot(timeSteps, statusUS);
title("Natick to GALILEO")
ylabel("# of satellites")
subplot(3,1,2);
plot(timeSteps, statusDE);
title("München to GALILEO")
ylabel("# of satellites")
subplot(3,1,3);
plot(timeSteps, statusIN);
title("Bangalore to GALILEO")
ylabel("# of satellites")


statusTable = [table(height(intervalsUS), height(intervalsDE), height(intervalsIN)); ...
    table(sum(intervalsUS.Duration)/3600, sum(intervalsDE.Duration)/3600, sum(intervalsIN.Duration)/3600); ...
    table(mean(intervalsUS.Duration/60), mean(intervalsDE.Duration/60), mean(intervalsIN.Duration/60)); ...
    table(mean(statusUS, 2), mean(statusDE, 2), mean(statusIN, 2)); ...
    table(min(statusUS), min(statusDE), min(statusIN)); ...
    table(max(statusUS), max(statusDE), max(statusIN))];
statusTable.Properties.VariableNames = ["Natick", "München", "Bangalore"];
statusTable.Properties.RowNames = ["Total # of intervals", "Total interval time (hrs)",...
    "Mean interval length (min)", "Mean # of satellites in view", ...
    "Min # of satellites in view", "Max # of satellites in view"];
statusTable



