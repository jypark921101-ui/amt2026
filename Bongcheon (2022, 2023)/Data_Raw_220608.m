clear;

filename_LICOR = ["Li-cor raw data\220608(1)_low conc..txt", "Li-cor raw data\220608(2)_high conc..txt"];

LICOR_Raw_1 = readtable(filename_LICOR(1));
LICOR_Raw_2 = readtable(filename_LICOR(2));

LICOR_Date_1 = LICOR_Raw_1.Date_Y_M_D_;
LICOR_TimeInt_1 = LICOR_Raw_1.Time_H_M_S_;
LICOR_Date_2 = LICOR_Raw_2.Date_Y_M_D_;
LICOR_TimeInt_2 = LICOR_Raw_2.Time_H_M_S_;

[YY, MM, DD] = ymd(LICOR_Date_1);
[HH, MI, SS] = hms(LICOR_TimeInt_1);
Time1 = datetime([YY, MM, DD, HH, MI, SS]);

[YY, MM, DD] = ymd(LICOR_Date_2);
[HH, MI, SS] = hms(LICOR_TimeInt_2);
Time2 = datetime([YY, MM, DD, HH, MI, SS]);

LICOR_Conc = [LICOR_Raw_1.CO2_ppm_; LICOR_Raw_2.CO2_ppm_];

LICOR = timetable([Time1; Time2], LICOR_Conc, 'VariableNames', "CO2");

LICOR_Ref_1 = [min(Time1), max(Time1)];
LICOR_Ref_2 = [min(Time2), max(Time2)];

Filenames = ["A", "B", "C", "D", "G", ...
    "H", "I", "J", "K",...
    "M", "N", "O", "P", "Q",...
    "R", "S", "T", "U", "X"];
FilePath = ["CO2 box raw data\1st calibration_low conc\", "CO2 box raw data\2nd calibration + 3rd measurement\"]; 
Box_Char = 'ABCDGHIJKMNOPQRSTUX'; 

LL = length(Filenames);

Box_1 = cell(LL, 1);
Box_2 = Box_1;

Box = {Box_1; Box_2};

Box_Ref_1 = NaT(LL, 2);
Box_Ref_2 = NaT(LL, 2);

for i = 1:LL
    Name = FilePath + Filenames(i) + ".xls";
    for j = 1:2
        if i == 19
            B1 = readtable(FilePath(j)+"X1.xls");
            B2 = readtable(FilePath(j)+"X2.xls");
            BB = [B1; B2];
        else
            BB = readtable(Name(j));
        end
        Box_Time_Raw = datetime(BB.Time, "InputFormat", 'MM-dd-yy HH:mm:ss');
        Box_Time_Start = Box_Time_Raw(1);
        Duration = seconds(Box_Time_Raw(end)-Box_Time_Raw(1));
        Time = dateshift(Box_Time_Start, "start", "second", (0:Duration)');
        Box_Conc = interp1(Box_Time_Raw, BB.CO2_PPM_, Time);
        Box_Temp = interp1(Box_Time_Raw, BB.Temp___, Time);
        Box_RH = interp1(Box_Time_Raw, BB.Humi__RH_, Time);
        Box{j}{i} = timetable(Time, Box_Conc, Box_Temp, Box_RH, 'VariableNames', ["CO2", "Temp", "RH"]);
        if j == 1
            Box_Ref_1(i, :) = [min(Time), max(Time)];
        else
            Box_Ref_2(i, :) = [min(Time), max(Time)];
        end
    end
end

save("220608_CO2_All_Raw", "Box", "Box_Char", "Box_Ref_1", "Box_Ref_2", "LICOR", "LICOR_Ref_1", "LICOR_Ref_2");