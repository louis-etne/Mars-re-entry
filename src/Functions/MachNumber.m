function mach = MachNumber(v, h, Atm)
    mach = v ./ SpeedOfSound(h, Atm); % No dimensions
end