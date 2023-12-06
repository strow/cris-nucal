% run_speccal_bygran.m
% ---------------------------------------------------------------------------------- 
% i = granule number index
% ---------------------------------------------------------------------------------- 
i = str2num(getenv('SLURM_ARRAY_TASK_ID'));

day = '210';
band = 1;
if band == 1;
   bandx = 'lw';
end
if band == 2;
   bandx = 'mw';
end
if band == 3;
   bandx = 'sw';
end

ifile = i;
sc = speccal_rtp_bygran(@sc_opts_for_speccal_bygran,i,band,day,[0 30]);
fout = ['/home/strow/Work/Cris/' day '/idps_day' day '_' bandx '_' int2str(ifile) ];
opts = sc_opts_for_speccal_bygran(band,day,[0 30],day);
save(fout,'opts','sc');



