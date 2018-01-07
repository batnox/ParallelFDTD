import os

compilation_string = "mpiicc -xMIC-AVX512 -qopenmp -fp-model fast=2 -ansi-alias xyz.c -lm -o exec/"
compilation_string_serial = "mpiicc -xMIC-AVX512 -fp-model fast=2 -ansi-alias xyz.c -lm -o exec/"

fin = open("fdtd_mpi_21.c")
lines =  fin.readlines()

for i in [256,512,1024,2048,4096]:
   lines[46] = "\tfptr=fopen(\"xstart_" + str(i) + ".dat\",\"r\");\n"

   fout = open("xyz.c","w")
   fout.writelines(lines)
   fout.close()

   os.system(compilation_string+"temp_exec/"+str(i)+"_mpi_openmp")
fin.close()

print "Hybrid Code Complied"
