% initKeplerUKF
function filter = initKeplerUKF(detection)

radarsphmeas = detection.Measurement;
[x, y, z] = sph2cart(deg2rad(radarsphmeas(1)),deg2rad(radarsphmeas(2)),radarsphmeas(3));
radarcartmeas = [x y z];
Recef2radar = detection.MeasurementParameters.Orientation;
ecefmeas = detection.MeasurementParameters.OriginPosition + radarcartmeas*Recef2radar;
% This is equivalent to:
% Ry90 = [0 0 -1 ; 0 1 0; 1 0 0]; % frame rotation of 90 deg around y axis
% nedmeas(:) = Ry90' * radarcartmeas(:);
% ecefmeas = lla2ecef(lla) +  nedmeas * dcmecef2ned(lla(1),lla(2));
initState = [ecefmeas(1); 0; ecefmeas(2); 0; ecefmeas(3); 0];

sigpos = 2;% m
sigvel = 0.5;% m/s^2

filter = trackingUKF(@keplerorbit,@cvmeas,initState,...
    'StateCovariance', diag(repmat([1000, 10000].^2,1,3)),...
    'ProcessNoise',diag(repmat([sigpos, sigvel].^2,1,3)));

end

function state = keplerorbit(state,dt)
% keplerorbit performs numerical integration to predict the state of
% Keplerian bodies. The state is [x;vx;y;vy;z;vz]

% Runge-Kutta 4th order integration method:
k1 = kepler(state);
k2 = kepler(state + dt*k1/2);
k3 = kepler(state + dt*k2/2);
k4 = kepler(state + dt*k3);

state = state + dt*(k1+2*k2+2*k3+k4)/6;

    function dstate=kepler(state)
        x =state(1,:);
        vx = state(2,:);
        y=state(3,:);
        vy = state(4,:);
        z=state(5,:);
        vz = state(6,:);

        mu = 398600.4405*1e9; % m^3 s^-2
        omega = 7.292115e-5; % rad/s
        
        r = norm([x y z]);
        g = mu/r^2;
        
        % Coordinates are in a non-intertial frame, account for Coriolis
        % and centripetal acceleration
        ax = -g*x/r + 2*omega*vy + omega^2*x;
        ay = -g*y/r - 2*omega*vx + omega^2*y;
        az = -g*z/r;
        dstate = [vx;ax;vy;ay;vz;az];
    end
end


% isDetectable 
function detectInput = isDetectable(tracker,time,covcon)

if ~isLocked(tracker)
    detectInput = zeros(0,1,'uint32');
    return
end
tracks = tracker.predictTracksToTime('all',time);
if isempty(tracks)
    detectInput = zeros(0,1,'uint32');
else
    alltrackid = [tracks.TrackID];
    isDetectable = zeros(numel(tracks),numel(covcon),'logical');
    for i = 1:numel(tracks)
        track = tracks(i);
        pos_scene = track.State([1 3 5]);
        for j=1:numel(covcon)
            config = covcon(j);
            % rotate position to sensor frame:
            d_scene = pos_scene(:) - config.Position(:);
            scene2sens = rotmat(config.Orientation,'frame');
            d_sens = scene2sens*d_scene(:);
            [az,el] = cart2sph(d_sens(1),d_sens(2),d_sens(3));
            if abs(rad2deg(az)) <= config.FieldOfView(1)/2 && abs(rad2deg(el)) < config.FieldOfView(2)/2
                isDetectable(i,j) = true;
            else
                isDetectable(i,j) = false;
            end
        end
    end
    
    detectInput = alltrackid(any(isDetectable,2))';
end
end


% assembleRadarInput 
function targets = assembleRadarInputs(station, constellation)
% For each satellite in the constellation, derive its pose with respect to
% the radar frame.
% inputs:
%          - station        :  LLA vector of the radar ground station
%          - constellation  :  Array of structs containing the ECEF position
%                              and ECEF velocity of each satellite
% outputs:
%          - targets        :  Array of structs containing the pose of each
%                              satellite with respect to the radar, expressed in the local
%                              ground radar frame (NED)

% Template structure for the outputs which contains all the field required
% by the radar step method
targetTemplate =  struct( ...
                'PlatformID', 0, ...
                'ClassID', 0, ...
                'Position', zeros(1,3), ...
                'Velocity', zeros(1,3), ...
                'Acceleration', zeros(1,3), ...
                'Orientation', quaternion(1,0,0,0), ...
                'AngularVelocity', zeros(1,3),...
                'Dimensions', struct( ...
                             'Length', 0, ...
                             'Width', 0, ...
                             'Height', 0, ...
                             'OriginOffset', [0 0 0]), ...
                'Signatures', {{rcsSignature}});


% First fill in the current satellite ECEF pose
targetPoses = repmat(targetTemplate, 1, numel(constellation));
for i=1:numel(constellation)
    targetPoses(i).Position = constellation(i).Position;
    targetPoses(i).Velocity = constellation(i).Velocity;
    targetPoses(i).PlatformID = constellation(i).PlatformID;
    % Orientation and angular velocity are left null, assuming satellite to
    % be point targets with a uniform rcs
end

% Then derive the radar pose in ECEF based on the ground station location
Recef2station = dcmecef2ned(station(1), station(2));
radarPose.Orientation = quaternion(Recef2station,'rotmat','frame');
radarPose.Position = lla2ecef(station);
radarPose.Velocity = zeros(1,3);
radarPose.AngularVelocity = zeros(1,3);

% Finally, take the difference and rotate each vector to the ground station
% NED axes
targets = targetPoses;
for i=1: numel(targetPoses)
    thisTgt = targetPoses(i);
    pos = Recef2station*(thisTgt.Position(:) - radarPose.Position(:));
    vel = Recef2station*(thisTgt.Velocity(:) - radarPose.Velocity(:)) - cross(radarPose.AngularVelocity(:),pos(:));
    angVel = thisTgt.AngularVelocity(:) - radarPose.AngularVelocity(:);
    orient = radarPose.Orientation' * thisTgt.Orientation;

    % Store into target structure array
    targets(i).Position(:) = pos;
    targets(i).Velocity(:) = vel;
    targets(i).AngularVelocity(:) = angVel;
    targets(i).Orientation = orient;
end
end


% addMeasurementParam
function dets = addMeasurementParams(dets, numdets, station)
% Add radar station information to the measurement parameters
Recef2station = dcmecef2ned(station(1), station(2));
for i=1:numdets
    dets{i}.MeasurementParameters.OriginPosition = lla2ecef(station);
    dets{i}.MeasurementParameters.IsParentToChild = true; % parent = ecef, child = radar
    dets{i}.MeasurementParameters.Orientation = dets{i}.MeasurementParameters.Orientation' * Recef2station;
end
end


% distanceFcn
function d = distanceFcn(track, truth)

estimate = track.State([1 3 5 2 4 6]);
true = [truth.Position(:) ; truth.Velocity(:)];
cov = track.StateCovariance([1 3 5 2 4 6], [1 3 5 2 4 6]);
d = (estimate - true)' / cov * (estimate - true);
end