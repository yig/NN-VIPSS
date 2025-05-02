
# replace the out_dir with your own output dir
# Figure 10 
./vipss -i ../../data/points/torus_n50.xyz -o {out_dir}/torus_n50.ply
./vipss -i ../../data/points/torus_halfsample.xyz -o {out_dir}/torus_halfsample.ply 
./vipss -i ../../data/points/bathtub_1k.xyz -o {out_dir}/bathtub_1k.ply 
./vipss -i ../../data/points/chair36.xyz -o {out_dir}/chair36.ply 

# Figure 11 
./vipss -i ../../data/points/torus_wires.xyz -o {out_dir}/torus_wires.ply
./vipss -i ../../data/points/doghead.xyz -o {out_dir}/doghead.ply
./vipss -i ../../data/points/hand_ok.xyz -o {out_dir}/hand_ok.ply
./vipss -i ../../data/points/walrus.xyz -o {out_dir}/walrus.ply -l 0.0005 --max_iter 1000

# Figure 11 
./vipss -i ../../data/points/Helmet.xyz -o {out_dir}/Helmet.ply
./vipss -i ../../data/points/brain_Ih.xyz -o {out_dir}/brain_Ih.ply
./vipss -i ../../data/points/Mobius.xyz -o {out_dir}/Mobius.ply

# Figure 13
./vipss -i ../../data/points/lord_quas.xyz -o {out_dir}/lord_quas.ply -l 0.0002 
./vipss -i ../../data/points/anchor.xyz -o {out_dir}/anchor.ply -l 0.001 

# Figure 1
./vipss -i ../../data/points/helmetMoustache.xyz -o {out_dir}/helmetMoustache.ply