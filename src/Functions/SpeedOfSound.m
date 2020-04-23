function speed = SpeedOfSound(T)
    addpath('Data');
    load('constants.mat', 'Atm');

    speed = sqrt(Atm.gamma * Atm.R * T); % m/s
end

