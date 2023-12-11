function sc = speccal_rtp(getopts,band,foutname,kmaxmin);
%---------------------------------------------------------------------------------
% Function speccal_rtp.m
%
% Read in a CrIS robs1 & rcalc RTP file and do spectral calibration.
%
% Calling Example:  speccal_rtp(@sc_fit_opts,3,filenameoutprefix,xtrack_max_min);
% 
% This cal says use params function "sc_fit_opts.m" and fit for band 3
%---------------------------------------------------------------------------------
% Main controlling options (for example, defined in sc_fit_opts.m)
% 
% induse:  channels to use for spectral calibration
% min_sc:  minimum BT contrast in selected channels
% max_dbt_surf:  maximum difference from simulated to observed surface BT
%---------------------------------------------------------------------------------

addpath ~/Git/matlib/h4tools          % rtpread
addpath ~/Git/matlib/cris             % cris_highres_chans
addpath /asl/matlab2012/cris/unapod  % apodization routines
addpath ~/Git/ccast/motmsc           % finterp2
addpath ~/Git/rtp_prod2/util         % rtp_sub_prof(prof,subset)
addpath ~/Git/matlib/aslutil          % radtobt
addpath ~/Git/matlib/time             % tai2dtime
%addpath /home/motteler/cris/cris_test/utils  % for now finterp2
addpath  ~/Git/ccast/motmsc/utils

% Limit and quantization of search
df = 0.625;
%dfmaxpercent =  2;
%dfsteppercent = 0.005;
dfmaxpercent =  8;
dfsteppercent = 0.020;

% Get proper frequencies for these data; Howard's routines
[n1,n2,n3,userLW,userMW,userSW, ichan] = cris_highres_chans(2);
f = cris_vchan(2, userLW, userMW, userSW);
nchan = n1+n2+n3;

% Subset to real channels
kstd = find(ichan <= nchan);
f = f(kstd);

lwi = 1:n1;
mwi = (n1+1):(n1+n2);
swi = (n1+n2+1):(n1+n2+n3);
% band indices used later
if band == 1;
   wi = lwi;
elseif band == 2
   wi = mwi;
elseif band == 3
   wi = swi;
end

% Fit params, changes with each run
% opts = sc_fit_opts(3);   % Save opts structure with output, flatten below
%                          % Note prof.solzen > 100 hard-coded below, fix later
opts = getopts(band,foutname,kmaxmin);   % Save opts structure with output, flatten below
                         % Note prof.solzen > 100 hard-coded below, fix later
induse       = opts.induse;
indsclo      = opts.indsclo;
indschi      = opts.indschi;
min_sc       = opts.min_sc;
max_dbt_surf = opts.max_dbt_surf;
indsurf      = opts.indsurf;
nmin         = opts.nmin;
glist        = opts.glist;
timeblock    = opts.timeblock;
relative     = opts.relative;
fout         = opts.fout;
szminmax     = opts.szminmax;

% Get info on function used in call to this function
opts.params = functions(getopts);

% Set up movement parameters for btcal
ninds = length(induse);
noffsets = 1 + 2*dfmaxpercent/dfsteppercent;
ffine = zeros(ninds,noffsets);
percentfine = zeros(1,noffsets);
for ii=1:noffsets
   percentfine(ii)= ((ii-1)*dfsteppercent - dfmaxpercent);
   d = percentfine(ii)/100;
   ffine(:,ii) = f(induse) + df*d;
end
clear d

% Obs counter set to 0
iblock = 0;
% Initialize mean obs counter   
for i=1:9; iblock(i) = 1; end
% Loop over input files in glist

for ifile=1:length(glist);

% Read current RTP file
   fin = glist(ifile).name;
   rtpdir = glist(ifile).folder;
   [head, hattr, prof, pattr] = rtpread(fullfile(rtpdir,fin));
   disp(['doing file ' fin])
   if isfield(prof,'rcalc');
      prof.rclr = prof.rcalc;
   end

% Reject profiles with any nan or obvious bad rcalc
   ibad = find(isnan(prof.rclr(1,:)) == 1 | prof.rclr(1,:) < 1E-5); 
   iok = setdiff(1:length(prof.nlevs),ibad);

   prof = rtp_sub_prof(prof,iok);
   clear iok ibad

% Sinc to Hamming if needed (or wanted)
   switch opts.robs_type
     case 'sinc'
       prof.robs1 = box_to_ham(prof.robs1);
   end
   switch opts.rcal_type
     case 'sinc'
%       prof.sarta_rclearcalc = box_to_ham(prof.rcalc);
      prof.rclr = box_to_ham(prof.rclr);
   end
% Get rid of guard channels so opts chan ids can be used directly     
   prof.robs1 = prof.robs1(kstd,:);
   prof.rclr = prof.rclr(kstd,:);

% Spectral contrast of observed spectra
   rlo = prof.robs1(indsclo,:);
   rhi = prof.robs1(indschi,:);
   dbto_sc = radtot(f(indschi),rhi) - radtot(f(indsclo),rlo); 

% Spectral contrast of calculated spectra
   rlo = prof.rclr(indsclo,:);
   rhi = prof.rclr(indschi,:);
   dbtc_sc = radtot(f(indschi),rhi) - radtot(f(indsclo),rlo); 

% Cloud effects for indsurf (if calcs are correct)
   bto_surf = radtot(f(indsurf),prof.robs1(indsurf,:));
   btc_surf = radtot(f(indsurf),prof.rclr(indsurf,:));

% Set up timeblocks
   mtime = tai2dtime(prof.rtime);
   tmin_obs = min(mtime);
   tmax_obs = max(mtime);
   
   if timeblock == 0
      d(1) = mtime(1);
      d(2) = mtime(end);
   else
      d = tmin_obs:minutes(timeblock):tmax_obs;
   end

% Loop over the d-1 time periods 
   for l = 2:length(d);
      xpercentfine = percentfine;
      for ifov = 1:9
         ip = find( mtime >= d(l-1) & mtime < d(l)           & ...
                    dbto_sc > min_sc & dbtc_sc > min_sc      & ...
                    abs(bto_surf - btc_surf) < max_dbt_surf  & ...
                    prof.ifov == ifov & ...
                    prof.solzen >= szminmax(1) & prof.solzen <= szminmax(2)   & ...
                    prof.xtrack >= opts.kmaxmin(1) & prof.xtrack < opts.kmaxmin(2) );
%                   prof.landfrac == 0 );
         ip5 = find( mtime >= d(l-1) & mtime < d(l)          & ...
                     dbto_sc > min_sc & dbtc_sc > min_sc     & ...
                     abs(bto_surf - btc_surf) < max_dbt_surf & ...
                    prof.ifov == 5 & ...
                    prof.solzen >= szminmax(1) & prof.solzen <= szminmax(2)   & ...
                    prof.xtrack >=opts.kmaxmin(1) & prof.xtrack <= opts.kmaxmin(2) );
         n = length(ip);

% Doppler correct band
         if opts.doppler
            psub = rtp_sub_prof(prof,ip);
            if length(ip) > 0
               dnu_ppm = doppler(psub);clear psub;
               sfac = 60;
               for i = 1:length(ip);
                  ro = prof.robs1(wi,ip(i));
                  [rad60 frq60] = finterp2(ro,f(wi),sfac);
                  frq60 = (1 + 1E-6*dnu_ppm(i))*frq60;
                  prof.robs1(wi,ip(i))= interp1(frq60(:),rad60,f(wi),'pchip','extrap');
               end
            end
         end

         
% Process current timeblock if sufficient FOVs
         if (n > nmin) 
            ro = nanmean(prof.robs1(induse,ip),2);
            if relative
               rc = nanmean(prof.robs1(induse,ip5),2);
            else
               rc = nanmean(prof.rclr(induse,ip),2);
            end
            bto = real(radtot(f(induse),ro));
            btc = real(radtot(f(induse),rc));
            sc(ifov).landfrac(iblock(ifov))     = mean(prof.landfrac(ip));
            sc(ifov).xtrack(iblock(ifov))       = mean(prof.xtrack(ip));
            sc(ifov).count(iblock(ifov))        = n;
            sc(ifov).tmin(iblock(ifov))         = d(l-1);
            sc(ifov).tmax(iblock(ifov))         = d(l);
            sc(ifov).rlat_mintime(iblock(ifov)) = min(prof.rlat(ip));
            sc(ifov).rlat_maxtime(iblock(ifov)) = max(prof.rlat(ip));
            sc(ifov).rlon_mintime(iblock(ifov)) = min(prof.rlon(ip));
            sc(ifov).rlon_maxtime(iblock(ifov)) = max(prof.rlon(ip));
            sc(ifov).bto(:,iblock(ifov))        = bto;
            sc(ifov).btc(:,iblock(ifov))        = btc;
            
% Here use finterp2 to create interpolated btc
% Then reinterp even finer below
            sfac = 60; % Way overkill
            [rad2, frq2] = finterp2(rc,f(induse),sfac);
            btc = rad2bt(frq2,rad2);
            
% Interpolate rcalc from f to ffine
            btcfine = interp1(frq2,btc,ffine,'pchip','extrap');

% Calculate correlation coefficients
            for ii = 1:noffsets
               junk = corrcoef(bto,btcfine(:,ii));
               junkall(ifile,ii) = junk(1,2);
               ccall(ii) = junk(1,2);
            end
            ii = find(ccall == max(ccall));
            if (length(ii) ~= 1)
               ii = ii(1);
               xpercentfine(ii(2:end)) = NaN;
            end
            sc(ifov).ccall(:,iblock(ifov))       = ccall;
            sc(ifov).corrcoef(iblock(ifov))      = ccall(ii);
            sc(ifov).percent_shift(iblock(ifov)) = xpercentfine(ii);
% Get std and mean of bias for best nu
            bias = bto - btcfine(:,ii);
            sc(ifov).btbias(iblock(ifov)) = nanmean(bias(:));
            sc(ifov).btbiasstd(iblock(ifov)) = nanstd(bias(:));
            iblock(ifov) = iblock(ifov)  + 1;
         end % (n > nmin)
      end % ifov
   end % timeblock loop
end % file loop

% Put results in matrix, using NaNs where no data

% First Get max lengths, and compute ppm shifts
mean_freq  = nanmean(f(induse));
channel_df = nanmean(diff(f(induse)));
for ifov = 1:9
   sc(ifov).ppm = 1E6 * sc(ifov).percent_shift * channel_df * 0.01 / mean_freq;
   sc_len(ifov) = length(sc(ifov).count);
end

ml = max(sc_len);
scall_ppm = nan(9,ml);
scall_corrcoef = nan(9,ml);
scall_count = nan(9,ml);
scall_mlat = nan(9,ml);
scall_xtrack = nan(9,ml);
scall_mtime = NaT(9,ml);
scall_mbto = nan(9,ml,length(induse));
scall_mbtc = nan(9,ml,length(induse));
scall_btbias = nan(9,ml);
scall_btbiasstd = nan(9,ml);

for ifov = 1:9
   l = 1:sc_len(ifov);
   scall_ppm(ifov,l) = sc(ifov).ppm;
   scall_corrcoef(ifov,l) = sc(ifov).corrcoef;
   scall_count(ifov,l) = sc(ifov).count;
   scall_xtrack(ifov,l) = sc(ifov).xtrack;
   scall_mlat(ifov,l) = 0.5 * (sc(ifov).rlat_mintime + sc(ifov).rlat_maxtime );
   scall_mlon(ifov,l) = 0.5 * (sc(ifov).rlon_mintime + sc(ifov).rlon_maxtime );
   scall_mtime(ifov,l) = sc(ifov).tmin + 0.5*minutes(opts.timeblock);
   scall_mbto(ifov,l,:) = sc(ifov).bto';
   scall_mbtc(ifov,l,:) = sc(ifov).btc';
   scall_btbias(ifov,l) =  sc(ifov).btbias;
   scall_btbiasstd(ifov,l) =  sc(ifov).btbiasstd;
end

save(fout,'opts', 'sc', 'scall_*'); 