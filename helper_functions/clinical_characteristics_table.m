clincal_charateristics_tables

 D = load ('alldat.mat');
 
 pivottable([D.SubjN, D.ID, D.gender, D.Age, D.handedness, D.paretic_side],...
     [],D.week,'length(unique(x))', 'datafilename','clinical_charateristics.txt');
 
 pivottable([D.control, D.week, D.ARAT_par, D.FM_total, D.FM_hand, ...
     D.MVC_FDI_par],[],D.week,'length(unique(x))', 'datafilename','clinical_scores.txt');