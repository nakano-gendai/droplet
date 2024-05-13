#======== make jobsqript & cmpile & conduct job ==========#
#!/bin/bash

#======== ジョブ変数入力 ==========#
#ジョブ名入力
echo "Input job name! (If blank, job_name = job)"
read JOB_NAME
if [ -z "$JOB_NAME" ]
then
	JOB_NAME='job'
fi
#ジョブ実行環境(small, small24VH, medium, large)
JOB_ENV='small'
#グループ番号
GROUP_NUM='23270'
#実行時間
ELA_TIME='00:15:00'
#VEノード数
VE_NODE='16'
#VEホスト数
VE_HOST='8'
#並列数
MPI_NUM='64'

#======== ジョブスクリプト作成 ==========#
touch ${JOB_NAME}.sh
echo "#!/bin/bash" > ${JOB_NAME}.sh
echo "#---- qsub option ----" >> ${JOB_NAME}.sh
echo "#PBS -q ${JOB_ENV}" >> ${JOB_NAME}.sh
echo "#PBS --group=${GROUP_NUM}" >> ${JOB_NAME}.sh
echo "#PBS -l elapstim_req=${ELA_TIME}" >> ${JOB_NAME}.sh
echo "#PBS --venode=${VE_NODE}" >> ${JOB_NAME}.sh
echo "#PBS --venum-lhost=${VE_HOST}" >> ${JOB_NAME}.sh
echo "#PBS -T necmpi" >> ${JOB_NAME}.sh
echo "#----- Program execution ----" >> ${JOB_NAME}.sh
echo "source /etc/profile.d/modules.sh" >> ${JOB_NAME}.sh
echo "module load NECSDK-sx/3.0.6" >> ${JOB_NAME}.sh
echo "module load 2decomp_fft-sx/1.5.847" >> ${JOB_NAME}.sh
echo "module load fftw-sx/3.3.8" >> ${JOB_NAME}.sh
echo "cd \$PBS_O_WORKDIR" >> ${JOB_NAME}.sh
echo "mpirun -np ${MPI_NUM} ./a.out" >> ${JOB_NAME}.sh

#=========　コンパイル & ジョブ投入　=============#
module load 2decomp_fft-sx/1.5.847
module load fftw-sx/3.3.8
module load NECSDK-sx/3.0.6
source /opt/nec/ve/bin/nlcvars.sh
mpinfort *.f90 -l2decomp_fft -lfftw3 -lasl_sequential -O3 && qsub ${JOB_NAME}.sh

if [ "$?" -eq "0" ]
then
	echo " "
	echo " ~~~~~  qsub succeeded!  ~~~~~"
	echo " "
	sstat
else
	echo " "
	echo 'qsub failed!'
	echo 'Check error message'
	echo " "
fi
