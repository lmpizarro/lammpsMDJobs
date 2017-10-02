atomsk --create fcc 4.5 U \
-duplicate 8 8 8 \
-select out sphere 0.5*box 0.5*box 0.5*box 10.0 \
-rmatom select \
esphere_in.cfg

atomsk --create fcc 4.5 Al \
-duplicate 8 8 8 \
-select in sphere 0.5*box 0.5*box 0.5*box 10.0 \
-rmatom select \
esphere_out.cfg


atomsk --merge 2 esphere_out.cfg esphere_in.cfg inclus_sphere.cfg
