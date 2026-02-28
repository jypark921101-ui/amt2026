clear;

LICOR_Raw = readtable("230208_LICOR.xlsx");
LICOR_Date = LICOR_Raw.Date_Y_M_D_;
LICOR_TimeInt = round(LICOR_Raw.Time_H_M_S_*86400);
[YY, MM, DD] = ymd(LICOR_Date);
HH = floor(LICOR_TimeInt/3600);
MI = floor(mod(LICOR_TimeInt, 3600)/60);
SS = mod(LICOR_TimeInt, 60);
Time = datetime([YY, MM, DD, HH, MI, SS]);
LICOR_Conc = LICOR_Raw.CO2_ppm_;

LICOR = timetable(Time, LICOR_Conc, 'VariableNames', "CO2");

Filenames = ["A", "C", "D", "F", "G", ...
    "H", "J", "K", "M", ...
    "N", "O", "P1", "Q1", "R1", ...
    "S1", "U1", "V1", "W1", "X1", "Z1"];
Box_Char = 'ACDFGHJKMNOPQRSUVWXZ';

LL = length(Filenames);

Box = cell(LL, 1);

for i = 1:LL
    BB = readtable(Filenames(i) + ".xls");
    Box_Time_Raw = datetime(BB.Time, "InputFormat", 'MM-dd-yy HH:mm:ss');
    Box_Time_Start = Box_Time_Raw(1);
    Duration = seconds(Box_Time_Raw(end)-Box_Time_Raw(1));
    Time = dateshift(Box_Time_Start, "start", "second", (0:Duration)');
    Box_Conc = interp1(Box_Time_Raw, BB.CO2_PPM_, Time);
    Box_Temp = interp1(Box_Time_Raw, BB.Temp___, Time);
    Box_RH = interp1(Box_Time_Raw, BB.Humi__RH_, Time);
    Box{i} = timetable(Time, Box_Conc, Box_Temp, Box_RH, 'VariableNames', ["CO2", "Temp", "RH"]);
end

save("230208_CO2_Raw", "Box", "Box_Char", "LICOR");