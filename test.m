day = '210';
band = '2';  % lw/mw/sw

i = 1;
sc = speccal_rtp_bygran(@sc_opts,band,'outputname_',[0 30]);
fout = ['./test_out/' day '/idps_day' day '_' bandx '_' int2str(ifile) ];
opts = sc_opts(band,'outputname',[0 30]);
save(fout,'opts','sc');



