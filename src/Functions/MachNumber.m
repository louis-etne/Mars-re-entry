function mach = MachNumber(v, h)
    mach = v ./ SpeedOfSound(h); % No dimensions
end