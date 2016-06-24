#PBS -N aSamplesFinal
#PBS -l nodes=1:ppn=4,walltime=96:00:00
#PBS -o /home/constantin/OEfiles
#PBS -e /home/constantin/OEfiles
#PBS -m abe
#PBS -M mailscluster@gmail.com
cd ~/matlab/projects/langevinCluster
/home/constantin/Software/Matlab/bin/matlab -nodesktop -nodisplay -nosplash -r "addLangevin ; clusterScript ; quit;"