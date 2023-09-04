import run_scf  # set the assignment


#run_scf.set_looper(100, 1.0e-10, 1.0e-8)
run_scf.set_looper_diis(100, 1.0e-10, 1.0e-8)