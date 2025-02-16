function [utc,jd] = Time(Period)
% UTC and juliandate calculation from period
% Input: Orbit's period
% Outputs:
% utc: Calculated UTC time from epoch (15-11-2023 12:00) to final time
% jd: Calculated juliandate from epoch (15-11-2023 12:00) to final time

t_0 = datetime('2023-11-15 12:00:00.000','InputFormat','yyyy-MM-dd HH:mm:ss.000'); % epoch is 15-11-2023 12:00:00
t_f = t_0 + seconds(Period-1); % Final time
time = t_0:seconds(1):t_f; % Creating array

jd = transpose(juliandate(time)); % Juliandate converter
%jd = jd(1:end-1);

% UTC time for calculations:
utc = datetime(jd,'ConvertFrom','juliandate'); utc.TimeZone = "UTC";
utc = datevec(utc);
end