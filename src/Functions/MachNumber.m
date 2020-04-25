function mach = MachNumber(v, T)
    mach = v ./ SpeedOfSound(T); % No dimensions
end