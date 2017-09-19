dials.import allow_multiple_sweeps=True th_8_1_000[1-3].cbf
dials.find_spots datablock.json nproc=3
dials.index datablock.json strong.pickle indexing.method=fft1d
dials.refine_bravais_settings indexed.pickle experiments.json 
dials.reindex indexed.pickle change_of_basis_op=a,b,c
dials.integrate bravais_setting_9.json reindexed_reflections.pickle nproc=3
