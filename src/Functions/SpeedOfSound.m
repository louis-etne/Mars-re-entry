function speed = SpeedOfSound(h)
    addpath('Data');
    load('constants.mat', 'Atm');

    speed = sqrt(Atm.gamma * Atm.R * Temperature(h)); % m/s
end

