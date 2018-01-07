import os

compilation_string = "icc -mmic -qopenmp -fp-model fast=2 -ansi-alias xyz.c -lm -o exec/"
compilation_string_serial = "icc -mmic -qopenmp -fp-model fast=2 -ansi-alias xyz.c -lm -o exec/"
'''
compilation_string = "icc -xMIC-AVX512 -qopenmp -fp-model fast=2 -ansi-alias xyz.c -lm -o exec/"
compilation_string_serial = "icc -xMIC-AVX512 -qopenmp -fp-model fast=2 -ansi-alias xyz.c -lm -o exec/"

compilation_string = "icc -qopenmp -fp-model fast=2 -ansi-alias xyz.c -lm -o exec/"
compilation_string_serial = "icc -qopenmp -fp-model fast=2 -ansi-alias xyz.c -lm -o exec/"

compilation_string = "gcc -qopenmp xyz.c -lm -o exec/"
compilation_string_serial = "gcc xyz.c -lm -o exec/"
'''

threads = [30, 32, 60, 64, 90, 120, 128, 160, 200, 240, 260, 272, 290]

fin = open("Block_Tiling/fdtd35_block_no_data_alignment.c")
lines =  fin.readlines()
fin1 = open("Block_Tiling/fdtd35_block_data_alignment_dynamic_array.c")
lines1 =  fin1.readlines()
fin2 = open("Block_Tiling/fdtd35_block_data_alignment_static_array.c")
lines2 =  fin2.readlines()

for i in [256,512,1024,2048,4096]:
    for j in [8,16,32,64]:
	lines[256] = "\t\tfptr=fopen(\"xstart_" + str(i) + ".dat\",\"r\");\n"
	lines[265] = "\t\titile=" + str(j) + ";\n"
	fout = open("xyz.c","w")
	fout.writelines(lines)
	fout.close()
	os.system(compilation_string+str(i)+"_guided_block_"+str(j)+"_no_data_alignment")
	
	lines1[256] = "\t\tfptr=fopen(\"xstart_" + str(i) + ".dat\",\"r\");\n"
	lines1[265] = "\t\titile=" + str(j) + ";\n"
	fout = open("xyz.c","w")
	fout.writelines(lines1)
	fout.close()
	os.system(compilation_string+str(i)+"_guided_dynamic_array_data_alignment_block_"+str(j))
	
	if i!=4096:
		lines2[50] = "\t\tfptr=fopen(\"xstart_" + str(i) + ".dat\",\"r\");\n"
		lines2[59] = "\t\titile=" + str(j) + ";\n"
		fout = open("xyz.c","w")
		fout.writelines(lines2)
		fout.close()
		os.system(compilation_string+str(i)+"_guided_static_array_data_alignment_block_"+str(j))

fin.close()
fin1.close()
fin2.close()

print "Block Code Complied"

fin = open("Row_Tiling/fdtd35_break_no_data_alignment.c")
lines =  fin.readlines()
fin1 = open("Row_Tiling/tiling_byte_counters.c")
lines1 = fin1.readlines()
fin2 = open("Row_Tiling/fdtd35_break_data_alignment_static_array.c")
lines2 =  fin2.readlines()
fin3 = open("Row_Tiling/fdtd35_break_data_alignment_dynamic_array.c")
lines3 =  fin3.readlines()

for i in [256,512,1024,2048,4096]:
	for j in [2,4,8,16]:
		for k in ["static", "dynamic", "guided"]:
			lines[254] = "\t\tfptr=fopen(\"xstart_" + str(i) + ".dat\",\"r\");\n"
			lines[650] = "\t\tint itile=" + str(j) + ";\n"
			lines[656] = "\t\t#pragma omp parallel for private(i,x,sine,sine1,j,omp2x,betax,const5x,const6x,extk,omp2y,betay,const6y,const5y,eytk) schedule(" + k + ")\n"
			fout = open("xyz.c","w")
			fout.writelines(lines)
			fout.close()
			os.system(compilation_string+str(i)+"_"+str(k)+"_tile_"+str(j)+"_no_data_alginment")
			
			if i!=4096:
				lines2[47] = "\t\tfptr=fopen(\"xstart_" + str(i) + ".dat\",\"r\");\n"
				lines2[441] = "\t\tint itile=" + str(j) + ";\n"
				lines2[447] = "\t\t#pragma omp parallel for private(i,x,sine,sine1,j,omp2x,betax,const5x,const6x,extk,omp2y,betay,const6y,const5y,eytk) schedule(" + k + ")\n"
				fout = open("xyz.c","w")
				fout.writelines(lines2)
				fout.close()
				os.system(compilation_string+str(i)+"_"+str(k)+"_static_array_data_alignment_tile_"+str(j))
		
			lines3[254] = "\t\tfptr=fopen(\"xstart_" + str(i) + ".dat\",\"r\");\n"
			lines3[650] = "\t\tint itile=" + str(j) + ";\n"
			lines3[656] = "\t\t#pragma omp parallel for private(i,x,sine,sine1,j,omp2x,betax,const5x,const6x,extk,omp2y,betay,const6y,const5y,eytk) schedule(" + k + ")\n"
			fout = open("xyz.c","w")
			fout.writelines(lines3)
			fout.close()
			os.system(compilation_string+str(i)+"_"+str(k)+"_dynamic_array_data_alignment_tile_"+str(j))
			
			lines1[256] = "\t\tfptr=fopen(\"xstart_" + str(i) + ".dat\",\"r\");\n"
			lines1[680] = "\t\tint itile=" + str(j) + ";\n"
			lines1[686] = "\t\t#pragma omp parallel for private(i,x,sine,sine1,j,omp2x,betax,const5x,const6x,extk,omp2y,betay,const6y,const5y,eytk) schedule(" + k + ")\n"
			fout = open("xyz.c","w")
			fout.writelines(lines1)
			fout.close()
			os.system(compilation_string+"/flops_"+str(i)+"_tile_"+str(j))

fin.close()
fin1.close()
fin2.close()

print "Loop Tiling Code Complied"
fin = open("Predefined_Tiling/fdtd20.c")
lines = fin.readlines()
fin1 = open("Predefined_Tiling/predefined_with_counters.c")
lines1 = fin1.readlines()
fin2 = open("Predefined_Tiling/fdtd20_static_array_data_alignment.c")
lines2 = fin2.readlines()

for i in [256,512,1024,2048,4096]:
	lines[248] = "\tfptr=fopen(\"xstart_"+str(i)+".dat\",\"r\");\n"
	lines[786] = "\t\t#pragma omp barrier\n"
	fout = open("xyz.c","w")
	fout.writelines(lines)
	fout.close()
	os.system(compilation_string+str(i)+"_barrier_predefined")
	
	lines1[249] = "\tfptr=fopen(\"xstart_"+str(i)+".dat\",\"r\");\n"
	fout = open("xyz.c","w")
	fout.writelines(lines1)
	fout.close()
	os.system(compilation_string+"/flops_"+str(i)+"_predefined")

	lines[786] = "\n"
	fout = open("xyz.c","w")
	fout.writelines(lines)
	fout.close()
	os.system(compilation_string+str(i)+"_no_barrier_predefined")
	
	lines2[42] = "\tfptr=fopen(\"xstart_"+str(i)+".dat\",\"r\");\n"
	lines2[576] = "\t\t#pragma omp barrier\n"
	fout = open("xyz.c","w")
	fout.writelines(lines2)
	fout.close()
	os.system(compilation_string+str(i)+"_barrier_predefined_static_array_data_alignment")

fin.close()
fin1.close()
fin2.close()

print "Predefined section Code Complied"
fin = open("Serial_Code/fdtd_serial.c")
lines = fin.readlines()
fin4 = open("Serial_Code/fdtd_serial_static_array_data_alignment.c")
lines4 = fin4.readlines()
fin1 = open("Naive_Code/naive_with_counters.c")
lines1 = fin1.readlines()
fin2 = open("Naive_Code/fdtd_naive.c")
lines2 = fin2.readlines()
fin3 = open("Naive_Code/fdtd_naive_static_data_alignment.c")
lines3 = fin3.readlines()

for i in [256,512,1024,2048,4096]:
	lines[255] = "\t\tfptr=fopen(\"xstart_" + str(i) + ".dat\",\"r\");\n"
	fout = open("xyz.c","w")
	fout.writelines(lines)
	fout.close()
	os.system(compilation_string_serial+str(i)+"_serial")
	
	if i!=4096:
		lines4[49] = "\t\tfptr=fopen(\"xstart_" + str(i) + ".dat\",\"r\");\n"
		fout = open("xyz.c","w")
		fout.writelines(lines4)
		fout.close()
		os.system(compilation_string_serial+str(i)+"_serial_static_array_data_alignment")

	lines1[316] = "\t\tfptr=fopen(\"xstart_" + str(i) + ".dat\",\"r\");\n"
	fout = open("xyz.c","w")
	fout.writelines(lines1)
	fout.close()
	os.system(compilation_string_serial+"/flops_"+str(i)+"_serial")
	
	lines2[235] = "\t\tfptr=fopen(\"xstart_" + str(i) + ".dat\",\"r\");\n"
	fout = open("xyz.c","w")
	fout.writelines(lines2)
	fout.close()
	os.system(compilation_string+str(i)+"_naive")
	
	lines3[29] = "\t\tfptr=fopen(\"xstart_" + str(i) + ".dat\",\"r\");\n"
	fout = open("xyz.c","w")
	fout.writelines(lines3)
	fout.close()
	os.system(compilation_string+str(i)+"_naive_static_array_data_alignment")

fin.close()
fin4.close()
fin1.close()
fin2.close()
fin3.close()
print "Serial Code Complied"
