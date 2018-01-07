import os

threads = [30, 32, 60, 64, 90, 120, 128, 160, 200, 240, 260, 272, 290]
data_size = [256,512,1024,2048,4096]
processes = [2,4,8,16,32]

for i in data_size:
	for j in processes:
		string = str(i)+"_mpi_"+str(j)
		for th in threads:
			os.system("mkdir exec/"+string+"_thread_"+str(th))
			os.system("cp exec/temp_exec/"+str(i)+"_mpi_openmp"+" exec/"+string+"_thread_"+str(th))
			os.system("cp exec/xstart_"+str(i)+".dat exec/"+string+"_thread_"+str(th))
			os.system("cp -r exec/anim_300 exec/"+string+"_thread_"+str(th))
			os.system("echo \"#PBS -N script\" > exec/"+string+"_thread_"+str(th)+"/jobscript")
			os.system("echo \"#PBS -l nodes="+str(j)+":knl7250:cache\n\" >> exec/"+string+"_thread_"+str(th)+"/jobscript")
			os.system("echo -e \'cd $PBS_O_WORKDIR\' >> exec/"+string+"_thread_"+str(th)+"/jobscript")
			os.system("echo \"export OMP_NUM_THREADS="+str(th)+"\" >> exec/"+string+"_thread_"+str(th)+"/jobscript")
			os.system("echo \"mpirun -machinefile \$PBS_NODEFILE -env OMP_NUM_THREADS "+str(th)+" ./"+str(i)+"_mpi_openmp"+" > "+string+"_threads_"+str(th)+"\" >> exec/"+string+"_thread_"+str(th)+"/jobscript")
			print "qsub -l walltime=00:60:00 -d ~/colfax/exec/"+string+"_thread_"+str(th)+" ~/colfax/exec/"+string+"_thread_"+str(th)+"/jobscript"

