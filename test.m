day = '210';
band = 2;  % lw/mw/sw

i = 1;
sc = speccal_rtp(@sc_opts,band,'outputname_',[0 30]);
fout = ['./test_out/output_band' num2str(band)];
opts = sc_opts(band,'outputname',[0 30]);
save(fout,'opts','sc');



