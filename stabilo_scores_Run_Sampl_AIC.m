function [output_bloc, v_devnorm] = stabilo_scores_JP(m_event_type_orig)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% scripte de JP (mla_grandstal_groupe2.m)
% v_devnorm is the normalized standard deviation (i.e. calculated on
% normalized RT, expressed in % of the median of correct trials) - the one
% we use to calculate PCT20 and PCT40
%
% INPUT :
%           m_event_type_orig : matrix, [code RT]
%                               Rt en ms
%
% MAJ :
% 30/01/18  RB
%       Correction de v_devnorm (stability sur fenetre glissante)
%       si 1 error ou miss apparait dans la fenetre => valeur par defaut =
%       100
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output_bloc   = {};


% les Rt doivent etre en ms
if mean(m_event_type_orig(:,2)) < 4
    m_event_type_orig(:,2) = m_event_type_orig(:,2)*1000;
end

m_event_type = m_event_type_orig; % this is the original one, before normalizing reaction times by the median, on non-zero values

output_bloc.nb_stim = size(m_event_type_orig,1);

% this is for normaization
v_rt_tool = m_event_type(:,2);

% keep only the good ones
v_x = m_event_type(:,1);
v_find_error = find(v_x >200);
v_find_ok    = find(v_x <200);
v_find_miss  = find(v_x >400);
v_find_false = find(v_x >200 & v_x <400);

output_bloc.pct_error = round(100*length(v_find_error)/length(v_x));
output_bloc.pct_ok    = round(100*length(v_find_ok)/length(v_x));
output_bloc.pct_miss  = round(100*length(v_find_miss)/length(v_x));
output_bloc.pct_false = round(100*length(v_find_false)/length(v_x));

output_bloc.avg_raw = nanmean(v_rt_tool);
output_bloc.med_raw = nanmedian(v_rt_tool);

output_bloc.avg_only_good = nanmean(v_rt_tool(v_find_ok));
output_bloc.med_only_good = nanmedian(v_rt_tool(v_find_ok));

s_averagetimepertrial_ms = nansum(v_rt_tool)+4000*length(v_find_error); % the sum of all reaction times, plus a penalty of 4s per error

output_bloc.avg_all_plus_penalization = s_averagetimepertrial_ms/length(v_x);



v_f_yesok = find(v_x==111);
v_rt_yesok = v_rt_tool(v_f_yesok);
s_med_yesok = round(median(v_rt_yesok));
v_f_nook = find(v_x==112);
v_rt_nook = v_rt_tool(v_f_nook);
s_med_nook = round(median(v_rt_nook));

v_f = find(v_x<200);
v_rt_tool = v_rt_tool(v_f);


v_ggg = find(v_rt_tool>0);
v_rt_tool = v_rt_tool(v_ggg);
s_med_rt = nanmedian(v_rt_tool); % we normalize by the median, KEEP THAT VALUE
s_med_rt_txt = round(s_med_rt);

m_event_type(:,2) = 100*m_event_type(:,2)/s_med_rt; %
% SO HERE, ALL THE REACTION TIMES HAVE BEEN NORMALIZED (DIVIDED BY THE
% MEDIAN). IT'S IN % OF THE MEDIAN
% SHALL WE DIVIDE BY THE MEAN INSTEAD ?
% QUESTION 1

%%
% we now compute our measure of stability
v_x = m_event_type(:,1);
v_rt = m_event_type(:,2);
% easy smoothing
v_ym = v_rt;
v_ys = zeros(size(v_ym));
for s_i = 2:(length(v_rt)-1)
%    if sum(v_x(s_i-1:s_i+1) > 200)<1
       v_ym(s_i) = mean(v_rt(s_i-1:s_i+1));
       v_ys(s_i) = std(v_rt(s_i-1:s_i+1),1);
%    else
%       v_ys(s_i) = 100;
%    end
end

% simply, each time there is a mistake, we set the std to 100, so that it
% cuts all the series
v_f = find(v_x>200);  % failures
v_ys(v_f) = 100; % set std to 100

v_ys(1) = 100;
v_ys(length(v_ys)) = 100;

v_devnorm = v_ys;
v_devnorm(1) = nan;
v_devnorm(end) = nan;



[output_bloc.pct20, output_bloc.pct40]     = Comput_PCT20_40(v_ys, m_event_type_orig);
[output_bloc.pct1500, output_bloc.pct3000] = Comput_PCT1500_3000(v_x, m_event_type_orig);



function [s_pct20, s_pct40] = Comput_PCT20_40(v_ys, m_event_type_orig)

v_barre = [1:40];
v_smalltick = (0:2:40);
v_bigtick = (0:10:40);

clear v_serix;

for s_i = 1:length(v_barre)
    v_test = (v_ys > v_barre(s_i)); % NOTE : v_t has already been modified in the last block, to be 3000 for each mistake
    v_pctmorestable(s_i) = 100*(length(v_ys)-sum(v_test))/(length(v_ys)-2); % pct of trials with a normalized std lower than the threshold (-2 because the first and the last points have erroneous std)
    v_f = find(v_test);
    v_diff = diff(v_f)-1; % that's the duration of all the series below v_barre(s_i)
    v_longeststable(s_i) = max(v_diff); % longest series with a normalized std lower than the threshold
    if (~isempty(find(v_diff)))
        v_tool = v_diff(find(v_diff))-4; % points start to accumulate after each series of 4 consecutive wins.  number 5 starts to bring points.
        v_tool = v_tool(find(v_tool>0));
        % v_tool contains all the points won for each good series
        v_serix(s_i) = sum(v_tool); % the number of points for that threshold
    else
        v_serix(s_i)=0;
    end;
end;

v_nbpoints_forpct40 = v_serix;

% the area under the curve AOC is
%v_serix = [0 v_serix(1:length(v_serix)-1)];

v_serix1 = v_serix(1:20); % we take only values less than 20%
s_aoc1 = sum(v_serix1);
s_maxaoc1 = (size(m_event_type_orig,1)-4)*length(v_serix1);
s_aocpct1 = 100*s_aoc1/s_maxaoc1;
s_pct20 = s_aocpct1;

s_aoc2 = sum(v_serix); % we take all
s_maxaoc2 = (size(m_event_type_orig,1)-4)*length(v_barre);  %
s_aocpct2 = 100*s_aoc2/s_maxaoc2;
s_pct40 = s_aocpct2;


function [s_pct1500, s_pct3000] = Comput_PCT1500_3000(v_x, m_event_type)

%m_event_type = m_event_type_orig; % this is the original one, before normalizing reaction times by the median, on non-zero values
% this is for normaization
v_rt = m_event_type(:,2);
v_f = find(v_x>200);  % failures
s_rttop = 4000;
v_rt(v_f) = s_rttop; % set to a reaction time of 4000
v_rt = [s_rttop; v_rt;s_rttop];

clear v_serix;

v_barre = [300:50:3000];

for s_i = 1:length(v_barre)
    v_test = (v_rt > v_barre(s_i)); % NOTE : v_t has already been modified in the last block, to be 3000 for each mistake
    v_f = find(v_test);
    v_diff = diff(v_f)-1; % that's the duration of all the series below v_barre(s_i)
    if (~isempty(find(v_diff)))
        v_tool = v_diff(find(v_diff))-4; % points start to accumulate after each series of 4 consecutive wins.  number 5 starts to bring points.
        v_tool = v_tool(find(v_tool>0));
        % v_tool contains all the points won for each good series
        v_serix(s_i) = sum(v_tool); % the number of points for that threshold
    else
        v_serix(s_i)=0;
    end;
end;

% the area under the curve AOC is
v_f = find(v_barre==1500); % we take only up to 1500 ms
v_serix1 = v_serix(1:v_f(1)); % we take only values less than 20%

s_aoc1 = sum(v_serix1);
s_maxaoc1 = (size(m_event_type,1)-4)*length(v_serix1);
s_pct1500 = 100*s_aoc1/s_maxaoc1;


s_aoc2 = sum(v_serix); % we take all
s_maxaoc2 = (size(m_event_type,1)-4)*length(v_barre);
s_pct3000 = 100*s_aoc2/s_maxaoc2;


