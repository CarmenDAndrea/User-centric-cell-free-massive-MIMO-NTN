function [theta_polar, phi_polar] = Geographic_to_polarConversion(elevation_angle, azimut_angle)
%elevation and azimuth angles in degrees

theta_polar=abs(elevation_angle)*pi/180;

phi_polar=azimut_angle*pi/180;

end