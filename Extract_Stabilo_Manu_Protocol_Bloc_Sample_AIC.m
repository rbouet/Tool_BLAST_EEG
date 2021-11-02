function output = Extract_Stabilo_Manu_Protocol_Bloc_Sample_AIC(path_file)



output = {};

% Load Manu's files
load(path_file)

output.trig_start_bloc = '';
output.Fs = HDR.SamplingFrequency;
output.ref_time_session = 0;   % en s
output.session_duration = HDR.NumberOfSamples/HDR.SamplingFrequency;
output.session_duration_without_eye = 0;
output.AIC.label = {};
output.AIC.code_label = [];
output.AIC.sample = [];
output.code_label = [];
output.sample = [];
output.OY.second = [];
output.OY.sample = [];
output.CY.second = [];
output.CY.sample = [];

output.bloc{1}.version = '+10ans';
output.bloc{1}.sample.Blast.lat = [1:5];    % en s
output.bloc{1}.sample.Blast.FB  = [71 71 71 71 71];
output.bloc{1}.sample.Blast.RT  = rand(1,5);    % en s
output.bloc{1}.sample.Blast.stab = [0 0 0 0 0];

output.bloc{1}.global.Blast.nb_stim = 5;
output.bloc{1}.global.Blast.pct_error = 0;
output.bloc{1}.global.Blast.pct_ok = 0;
output.bloc{1}.global.Blast.pct_miss = 0;
output.bloc{1}.global.Blast.pct_false = 0;
output.bloc{1}.global.Blast.avg_raw = 0;
output.bloc{1}.global.Blast.med_raw = 0;
output.bloc{1}.global.Blast.avg_only_good = 0;
output.bloc{1}.global.Blast.med_only_good = 0;
output.bloc{1}.global.Blast.avg_all_plus_penalization = 0;
output.bloc{1}.global.Blast.pct20 = 0;
output.bloc{1}.global.Blast.pct40 = 0;
output.bloc{1}.global.Blast.pct1500 = 0;
output.bloc{1}.global.Blast.pct3000 = 0;

output.bloc{1}.duration = output.session_duration;   % de la première stim à la dernière réponse
output.bloc{1}.start_stop = [0 output.session_duration]; 
output.bloc{1}.global.AIC.aic_general_nb = 0;
output.bloc{1}.global.AIC.aic_general_pct = 0;
output.bloc{1}.global.AIC.aic_focal_nb = 0;
output.bloc{1}.global.AIC.aic_focal_pct = 0;
output.bloc{1}.global.AIC.other_nb = 0;
output.bloc{1}.global.AIC.other_pct = 0;
output.bloc{1}.global.AIC.seiz_general_nb = 0;
output.bloc{1}.global.AIC.seiz_general_pct = 0;
output.bloc{1}.global.AIC.seiz_focal_nb = 0;
output.bloc{1}.global.AIC.seiz_focal_pct = 0;
    
output.bloc{1}.sample.AIC.befor_stim.nb_aic = [0 0 0 0 0];
output.bloc{1}.sample.AIC.befor_stim.cumultime_aic = [0 0 0 0 0];
output.bloc{1}.sample.AIC.befor_stim.cumulPCT_aic = [0 0 0 0 0];
output.bloc{1}.sample.AIC.nb_aic = [0 0 0 0 0];
output.bloc{1}.sample.AIC.aic_duration = [0 0 0 0 0];
output.bloc{1}.sample.AIC.aic_pct_duration = [0 0 0 0 0];
output.bloc{1}.sample.AIC.aic_position_stim = {[],[],[],[],[]};
output.bloc{1}.sample.AIC.type_aic = {[],[],[],[],[]};
     
output.bloc{1}.focus.m10.serie_max = 1;
output.bloc{1}.focus.m10.lat = 1;
output.bloc{1}.focus.m10.serie_dur = 1;
output.bloc{1}.focus.m10.total_dur = 1;
output.bloc{1}.focus.m20.serie_max = 1;
output.bloc{1}.focus.m20.lat = 2;
output.bloc{1}.focus.m20.serie_dur = 1;
output.bloc{1}.focus.m20.total_dur = 1;
output.bloc{1}.focus.m100.serie_max = 1;
output.bloc{1}.focus.m100.lat = 3;
output.bloc{1}.focus.m100.serie_dur = 1;
output.bloc{1}.focus.m100.total_dur = 1;
   
