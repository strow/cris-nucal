day = '210';
band = mw;  % lw/mw/sw

sc = speccal_rtp_bygran(@sc_opts,i,band,day,[0 30]);
fout = ['./test_out/' day '/idps_day' day '_' bandx '_' int2str(ifile) ];
opts = sc_opts(band,day,[0 30],day);
save(fout,'opts','sc');



