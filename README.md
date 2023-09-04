# directjk
linked Libcint for int1e & int2e  
write by C and python  
谨以此纪念蔚欣陪我度过的快乐时光  
祝愿她美国留学生活顺利、愉快，身体健康  
希望她两年后Python写的比我好！  
gcc -shared -o libcalc_int2e.so -fPIC calculate_int2e.c -O3 -lcint   
gcc -shared -o libcalc_int1e.so -fPIC calculate_int1e.c -O3 -lcint  
python setup.py build_ext --inplace
