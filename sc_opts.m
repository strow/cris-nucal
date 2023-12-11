function opts = sc_opts(band,run_name_prefix,kmaxmin);
%---------------------------------------------------------------------------------
%run_name_prefix = 'mw_clear_'';
%---------------------------------------------------------------------------------
% Likely global params good for all three bands

%datadir = '/asl/rtp/rtp_cris2_ccast_hires_noaa_a4/clear/2018';
%datadir = '/asl/rtp/rtp_cris2_ccast_hires_a2v4_ref/clear/2018';
%datadir = '/asl/rtp/rtp_cris2_ccast_hires_j1v4_a2v4/clear/2018';
%  datadir = '/asl/rtp/rtp_cris_NOAA_sdr_lowres/allfov/2019/060';


%datadir = '/asl/rtp/cris/snpp_sideswitch_sdr_hires/allfov/2019/060'
%datadir = '/asl/rtp/cris/snpp_sideswitch_sdr_hires/allfov/2021/195';
datadir = '/asl/rtp/cris/npp_ADL_hires/allfov/2021/219';

datadir = '/asl/rtp/cris/snpp_072023restart/allfov/2023/255';
datadir = '/asl/rtp/cris/j01_noaa/allfov/2023/275';;
datadir = 'Data';
%datadir = '/asl/rtp/cris/j02_calval/allfov/2023/041';
%datadir = '/asl/rtp/rtp_cris_ccast_hires/allfov/2018/061';
%datadir = '/asl/rtp/rtp_cris2_ccast_hires_def_j1m1/clear/2018/';
%datadir = '/asl/rtp/rtp_cris2_ADL_jan21_sdr_hires/clear/2018/008/';
%datadir = '/asl/rtp/rtp_cris2_ccast_hires_j1v3_a2v3/clear/2018/';
% datadir = '/asl/rtp/rtp_cris2_ccast_hires_j1v3_UWa2/clear/2018/'
% clist = { ...
%  'cris_ecmwf_csarta_clear_d20180121.rtp' };
% 'cris_ecmwf_csarta_clear_d20180107.rtp'  ... 
% 'cris_ecmwf_csarta_clear_d20180107.rtp'  ... 
% 'cris_ecmwf_csarta_clear_d20180108.rtp'  ... 
% 'cris_ecmwf_csarta_clear_d20180109.rtp'  ... 
% 'cris_ecmwf_csarta_clear_d20180110.rtp'  ... 
% 'cris_ecmwf_csarta_clear_d20180111.rtp'  ... 
% 'cris_ecmwf_csarta_clear_d20180112.rtp'  ... 
% 'cris_ecmwf_csarta_clear_d20180113.rtp'  ... 
% 'cris_ecmwf_csarta_clear_d20180114.rtp'  ... 
% 'cris_ecmwf_csarta_clear_d20180115.rtp'  ... 
% 'cris_ecmwf_csarta_clear_d20180116.rtp'};
% 
% for i=1:length(clist)
%    glist(i).name   = clist{i};
%    glist(i).folder = datadir;
% end

% Recursive list
%glist = dir([datadir '**/cris_ecmwf_csarta_clear*.rtp']);
%glist = dir(fullfile(datadir,'cris2_ecmwf_csarta_clear_d201803*.rtp'));
%glist = dir(fullfile(datadir,'cris2_ecmwf_csarta_clear_d2018*.rtp'));
%glist = dir(fullfile(datadir,'cris_ecmwf_csarta_allfov_d*.rtp'));

%glist = dir(fullfile(datadir,'cris_sdr_ecmwf_csarta_allfov_d*.rtp'));
glist = dir(fullfile(datadir,'cris_ecmwf_csarta_allfov_d*.rtp'));
%glist = dir(fullfile(datadir,'cris_sdr_ecmwf_csarta_clear_d*.rtp'));

opts.glist = glist;

% % These values appropriate for csarta files
opts.robs_type = 'sinc';
opts.rcal_type = 'notsinc';

% % These values appropriate for isarta files
% opts.robs_type = 'sinc';
% opts.rcal_type = 'sinc';

opts.indsurf = 501;  % Check for clouds/errors using 964 cm-1
opts.nmin = 3;    %WAS 10  % Min number of obs in a time period
opts.timeblock = 0;  % Average 6 minutes of data at a time

opts.relative = false;  % Correlation vs calc, if true, correlation vs FOV5
opts.kmaxmin = kmaxmin;
opts.doppler = false;
%opts.szminmax = [90 180];
%---------------------------------------------------------------------------------
if band == 1
% Strat CO2 (not used)
% induse = [3:20,34:72]; % Strat, avoiding Q-branches
% indsclo = 54; % 683.1 wn @ 218K
% indschi = 53; % 682.5 wn @ 238K
% min_sc = 10;
% max_dbt_surf = 10;

% Band1 mid/upper trop CO2
%  induse = [100:162]; % mid/upper trop CO2
%  indsclo = 124; % 726.3 wn & 265K
%  indschi = 123; % 726.9 wn @ 235K
%  min_sc = 5;%10;
%  max_dbt_surf = 5;

% Water lines in window region, used due to stable background emission
   opts.induse = 200:280;
   opts.indsclo = 239; % 726.3 wn & 265K  window
   opts.indschi = 240; % 726.9 wn @ 235K  window
   opts.min_sc  = 1;  %5%10;
   opts.max_dbt_surf = 5;
   opts.fout = [run_name_prefix '_lw'];
%   opts.szminmax = [1 180];
   opts.szminmax = [90 180];  % night
%   opts.szminmax = [0 180];
end
%---------------------------------------------------------------------------------
if band == 2
% Band2 mid/upper trop water
   opts.induse = [1164:1307];
   opts.indsclo = 1188; 
   opts.indschi = 1182; 
   opts.min_sc = 3;
   opts.max_dbt_surf = 5;
   opts.szminmax = [90 180];  % night
%   opts.szminmax = [0 180];
   opts.fout = [run_name_prefix '_mw'];
end
%---------------------------------------------------------------------------------
if band == 3
% 2155-2370 wn upper strat
%    opts.induse=1856:1883; 
%    opts.indschi = 1863;
%    opts.indsclo = 1864;
%    opts.min_sc = 0.5;
%    opts.max_dbt_surf = 20;  % Usually 5, but for these data just doesn't matter
%    opts.fout = [run_name_prefix 'strat_sw'];
%end
% % % % CO region
   opts.induse= 1593:1618;
   opts.indschi = 1600;
   opts.indsclo = 1602;
   opts.min_sc = 0.5;
   opts.szminmax = [90 180];  % night
   opts.max_dbt_surf = 5;  % Usually 5, but for these data just doesn't matter
   opts.fout = [run_name_prefix 'co_sw'];
% end
% All nu region
%    opts.induse= 1579:2211;
%    opts.indschi = 1600;
%    opts.indsclo = 1602;
%    opts.min_sc = 0.1;
%    opts.max_dbt_surf = 5;  % Usually 5, but for these data just doesn't matter
%    opts.fout = [run_name_prefix 'sw_allchans'];
end
%---------------------------------------------------------------------------------
