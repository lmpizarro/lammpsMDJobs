atomsk --create fcc 4.5 U \
-duplicate 10 10 10 \
-select out sphere 0.5*box 0.5*box 0.5*box 12.0 \
-rmatom select \
esphere_in.cfg

atomsk --create fcc 4.5 Al \
-duplicate 10 10 10 \
-select in sphere 0.5*box 0.5*box 0.5*box 12.0 \
-rmatom select \
esphere_out.cfg


atomsk --merge 2 esphere_out.cfg esphere_in.cfg sphere_in.cfg


rm esphere_in.cfg esphere_out.cfg


atomsk sphere_in.cfg sphere_in.lmp

