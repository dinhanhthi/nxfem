function changeResolution(dim)
    s = settings;s.matlab.desktop.DisplayScaleFactor
    s.matlab.desktop.DisplayScaleFactor.PersonalValue = dim;
end