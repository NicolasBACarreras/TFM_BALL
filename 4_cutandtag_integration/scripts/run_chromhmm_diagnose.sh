module load java

samples_directory=/
chr_length_file=/gpfs/projects/bsc08/shared_projects/IJC_3Dchromatin/genomes/human/GRCh38/current/Homo_sapiens.GRCh38.104.dna.all_chr.fa.fai
cell_mark_file=/home/bsc59/bsc59153/chromhmm_patient2/cell_mark_diagnose.txt
output_binary_directory=/home/bsc59/bsc59153/chromhmm_patient2/BALL_binary_directory_diagnose
output_model_directory=/home/bsc59/bsc59153/chromhmm_patient2/BALL_diagnose_model
num_states=13
java -mx4000M -jar /gpfs/projects/bsc08/shared_projects/IJC_3Dchromatin/programs/ChromHMM_1_24/ChromHMM.jar BinarizeBam -paired $chr_length_file $samples_directory $cell_mark_file $output_binary_directory
java -mx20000M -jar /gpfs/projects/bsc08/shared_projects/IJC_3Dchromatin/programs/ChromHMM_1_24/ChromHMM.jar LearnModel -p 48 $output_binary_directory $output_model_directory $num_states GRCh38
