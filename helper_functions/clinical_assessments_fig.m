clinical_assessment_fig

D = load ('alldat.mat');

D.mvc_ratio = D.MVC_FDI_unpar ./ D.MVC_FDI_par

figure;
lineplot(D.week,D.ARAT_par,'split', D.SubjN,'style_symbols4*2');
xlabel ('Week')
ylabel ('ARAT lesioned side')
title ('Patients')

figure;
lineplot(D.week,D.FM_total,'split', D.SubjN,'style_symbols4*2');
xlabel ('Week')
ylabel ('FM_total lesioned side')
title ('Patients')

figure;
lineplot(D.week,D.FM_hand,'split', D.SubjN,'style_symbols4*2');
xlabel ('Week')
ylabel ('FM_hand lesioned side')
title ('Patients')

figure;
lineplot(D.week,D.mvc_ratio,'split', D.SubjN,'style_symbols4*2');
xlabel ('Week')
ylabel ('FM_hand lesioned side')
title ('Patients')

