function e = Oswald(AR, sweep)
if sweep < deg2rad(30)
    e = 1.78*(1 - 0.044*AR.^0.68) - 0.64;
else
    e = 4.61*(1 - 0.044*AR.^0.68).*cos(sweep).^0.15 - 3.1;
end
end