module load java

samples_directory=/
chr_length_file=/gpfs/projects/bsc08/shared_projects/IJC_3Dchromatin/genomes/human/GRCh38/current/Homo_sapiens.GRCh38.104.dna.all_chr.fa.fai
cell_mark_file=/home/bsc59/bsc59153/chromhmm_patient2/cell_mark_relapse.txt
output_binary_directory=/home/bsc59/bsc59153/chromhmm_patient2/BALL_binary_directory_relapse_16
output_model_directory=/home/bsc59/bsc59153/chromhmm_patient2/BALL_relapse_model_16
num_states=16
java -mx4000M -jar /gpfs/projects/bsc08/shared_projects/IJC_3Dchromatin/programs/ChromHMM_1_24/ChromHMM.jar BinarizeBam -paired $chr_length_file $samples_directory $cell_mark_file $output_binary_directory
java -mx20000M -jar /gpfs/projects/bsc08/shared_projects/IJC_3Dchromatin/programs/ChromHMM_1_24/ChromHMM.jar LearnModel -p 48 $output_binary_directory $output_model_directory $num_states GRCh38
