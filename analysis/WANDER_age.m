slist = [1:5 8:13 15:20 22:26];

% year,month,day; sex: 0 = female
DOB{1} = datenum(1992,04,28);  
MEG{1} = datenum(2016,05,09);
SEX(1) = 1;
% year,month,day
DOB{2} = datenum(1991,04,09);  
MEG{2} = datenum(2016,05,10);
SEX(2) = 0;
% year,month,day
DOB{3} = datenum(1986,02,17);  
MEG{3} = datenum(2016,05,10);
SEX(3) = 0;
% year,month,day
DOB{4} = datenum(1992,06,08);  
MEG{4} = datenum(2016,05,11);
SEX(4) = 1;
% year,month,day
DOB{5} = datenum(1989,02,17);  
MEG{5} = datenum(2016,05,11);
SEX(5) = 0;
% year,month,day
DOB{6} = datenum(1993,06,30);  
MEG{6} = datenum(2016,05,12);
SEX(6) = 1;
% year,month,day
DOB{7} = datenum(1992,02,24);  
MEG{7} = datenum(2016,05,13);
SEX(7) = 0;
% year,month,day
DOB{8} = datenum(1989,09,25);  
MEG{8} = datenum(2016,05,13);
SEX(8) = 0;
% year,month,day
DOB{9} = datenum(1989,04,21);  
MEG{9} = datenum(2016,05,18);
SEX(9) = 0;
% year,month,day
DOB{10} = datenum(1994,07,10);  
MEG{10} = datenum(2016,05,20);
SEX(10) = 1;
% year,month,day
DOB{11} = datenum(1994,11,28);  
MEG{11} = datenum(2016,05,23);
SEX(11) = 1;
% year,month,day
DOB{12} = datenum(1994,06,10);  
MEG{12} = datenum(2016,05,23);
SEX(12) = 0;
% year,month,day
DOB{13} = datenum(1993,12,16);  
MEG{13} = datenum(2016,05,24);
SEX(13) = 1;
% year,month,day
DOB{14} = datenum(1990,03,30);  
MEG{14} = datenum(2016,05,30);
SEX(14) = 1;
% year,month,day
DOB{15} = datenum(1990,10,04);  
MEG{15} = datenum(2016,06,21);
SEX(15) = 1;
% year,month,day
DOB{16} = datenum(1989,12,08);  
MEG{16} = datenum(2016,06,02);
SEX(16) = 0;
% year,month,day
DOB{17} = datenum(1996,01,06);  
MEG{17} = datenum(2016,06,03);
SEX(17) = 1;
% year,month,day
DOB{18} = datenum(1994,09,19);  
MEG{18} = datenum(2016,06,03);
SEX(18) = 0;
% year,month,day
DOB{19} = datenum(1995,02,24);  
MEG{19} = datenum(2016,06,06);
SEX(19) = 0;
% year,month,day
DOB{20} = datenum(1994,09,02);  
MEG{20} = datenum(2016,06,08);
SEX(20) = 1;
% year,month,day
DOB{21} = datenum(1990,10,21);  
MEG{21} = datenum(2016,06,10);
SEX(21) = 0;
% year,month,day
DOB{22} = datenum(1990,10,26);  
MEG{22} = datenum(2016,06,10);
SEX(22) = 0;
% year,month,day
DOB{23} = datenum(1990,11,03);  
MEG{23} = datenum(2016,06,13);
SEX(23) = 1;
% year,month,day
DOB{24} = datenum(1993,05,02);  
MEG{24} = datenum(2016,06,14);
SEX(24) = 0;
% year,month,day
DOB{25} = datenum(1991,06,06);  
MEG{25} = datenum(2016,06,17);
SEX(25) = 1;
% year,month,day
DOB{26} = datenum(1994,05,27);  
MEG{26} = datenum(2016,06,22);
SEX(26) = 1;

for isubject = slist
    age(isubject) = abs(MEG{isubject} - DOB{isubject}) / 365;
end
age_avg = mean(age(slist));
age_std = std(age(slist));
age_min = min(age(slist));
age_max = max(age(slist));
nr_male = sum(SEX(slist));
nr_female = length(slist) - sum(SEX(slist));
