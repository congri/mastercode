#PBS -N temperedOptim
#PBS -l nodes=1:ppn=2,walltime=120:00:00
#PBS -o /home/constantin/OEfiles
#PBS -e /home/constantin/OEfiles
#PBS -m abe
#PBS -M mailscluster@gmail.com
cd ~/matlab/projects/mastercode
/home/constantin/Software/Matlab/bin/matlab -nodesktop -nodisplay -nosplash -r "adddir ; main ; quit;"