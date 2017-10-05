t=0 ; 
for s in dump.*.jpg; do 
mv $s dump$t.jpg ; t=`expr $t + 1`; 
done

