function speed = SpeedOfSound(h, Atm)
    speed = sqrt(Atm.gamma * Atm.R * Temperature(h)); % m/s
end

