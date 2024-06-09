import matplotlib.pyplot as plt

# Read in probability file and set threshold 

threshold_high = 0.75
threshold_low = 0.2
state_dictionary = {}
translator = {}
reso = 20000

with open('../data/Reference_matrix_Ub.csv', 'r') as fh:
    next(fh)
    for line in fh:
        fields = line.strip().split(',')
        values = fields[0]
        key = ','.join(fields[1:7])  # Joining the six fields with commas
        state_dictionary[key] = values
        #print("Key:", key)
        #print("Values:", values)
        
        
with open('../chromhmm_outputs/chromhmm_15_total_states/BALL_total_model/emissions_15.txt', 'r') as fh:
    next(fh)  # This skips the first line
    for line in fh:
    
        i, mark1, mark2, mark3, mark4, mark6, mark5 = line.split("\t")
        
        if float(mark1) < 0.07:
            mark1 = 0
            
        if float(mark2) < 0.07:
            mark2 = 0
            
        if float(mark3) < 0.07:
            mark3 = 0
            
        if float(mark4) < 0.07:
            mark4 = 0
            
        if float(mark5) < 0.07:
            mark5 = 0
        
        if float(mark6) < 0.07:
            mark6 = 0
            
           
        max_value = max(float(mark1), float(mark2), float(mark3), float(mark4), float(mark5), float(mark6))
        #print(max_value)
        
        
        """         if i == "4":
            print(mark1, mark2, mark3, mark4, mark5, mark6)
        if i == "5":
            print(mark1, mark2, mark3, mark4, mark5, mark6) """
        
        
        if max_value == 0:
            
            pattern = "F,F,F,F,F,F"
            
            try:    
                print(f"E{i}: {state_dictionary[pattern]}")
                translator[f"E{i}"] = state_dictionary[pattern]
            except KeyError:
                print(f"E{i} has no translation: {pattern}")
                
            continue

 
        #print(float(mark1)/float(max_value), float(mark2)/float(max_value), float(mark3)/float(max_value), float(mark4)/float(max_value), float(mark5)/float(max_value), float(mark6)/float(max_value))

        

        


        if  float(mark1)/float(max_value) > threshold_low:
            mark1 = "T"
        else:
            mark1 = "F"
        
        if  float(mark2)/float(max_value) > threshold_low:
            mark2 = "T"
        else:
            mark2 = "F"        
        
        if  float(mark3)/float(max_value) > threshold_low:
            mark3 = "T"
        else:
            mark3 = "F"        
        
        if  float(mark4)/float(max_value) > threshold_low:
            mark4 = "T"
        else:
            mark4 = "F"        
        
        if  float(mark5)/float(max_value) > threshold_low:
            mark5 = "T"
        else:
            mark5 = "F"
            
        if float(mark6)/float(max_value) > threshold_low:
            mark6 = "T"
        else:
            mark6 = "F"
        
        #print(f"{mark1},{mark2},{mark3},{mark4},{mark5},{mark6}")
        #pattern = f"{mark6},{mark5},{mark4},{mark3},{mark1},{mark2}"
        pattern = f"{mark1},{mark2},{mark3},{mark4},{mark5},{mark6}"   #This for 15 states total model
        #print(pattern)

        try:    
            print(f"E{i}: {state_dictionary[pattern]} -> {pattern}")
            translator[f"E{i}"] = state_dictionary[pattern]
        except KeyError:
            print(f"E{i} has no translation: {pattern}. Find manual annotation")

##########################################################################################################################

#Manually annotate missing states

translator["E4"] = "EPr"
""" translator[f"E10"] = "PoisedP"
translator[f"E11"] = "PolycombRepressed"
translator[f"E6"] = "EA"
translator[f"E7"] = "EA"
translator[f"E8"] = "LowSignal" """



########################################################################################################################## 

# Initialize dictionaries to store state counts
state_counts_relapse = {}
translated_state_counts_relapse = {}

with open('../chromhmm_outputs/chromhmm_15_total_states/BALL_total_model/relapse_15_segments.bed', 'r') as fh_in, \
     open('../results/relapse_15_segments_translated.bed', 'w') as fh_out:

    for line in fh_in:
        chrom, pos1, pos2, state = line.split()
        
        # Translate state
        state_translated = translator.get(state, state)  # Use translator.get() to handle cases where state is not found
        
        # Update counts for translated and original states
        try:
            translated_state_counts_relapse[state_translated] += float(int(pos2) - int(pos1))
        except KeyError:
            translated_state_counts_relapse[state_translated] = float(int(pos2) - int(pos1))
            
        try:
            state_counts_relapse[state] += float(int(pos2) - int(pos1))
        except KeyError:
            state_counts_relapse[state] = float(int(pos2) - int(pos1))
        
        # Write to the new BED file with translated state
        fh_out.write(f"{chrom}\t{pos1}\t{pos2}\t{state_translated}\n")
        
        
        

        
        #print(chrom, pos1, pos2, state, state_translated)
    
# Initialize dictionaries to store state counts
state_counts_diagnose = {}
translated_state_counts_diagnose = {}

with open('../chromhmm_outputs/chromhmm_15_total_states/BALL_total_model/diagnose_15_segments.bed', 'r') as fh_in, \
     open('../results/diagnose_15_segments_translated.bed', 'w') as fh_out:

    for line in fh_in:
        chrom, pos1, pos2, state = line.split()
        
        # Translate state
        state_translated = translator.get(state, state)  # Use translator.get() to handle cases where state is not found
        
        # Update counts for translated and original states
        try:
            translated_state_counts_diagnose[state_translated] += float(int(pos2) - int(pos1))
        except KeyError:
            translated_state_counts_diagnose[state_translated] = float(int(pos2) - int(pos1))
            
        try:
            state_counts_diagnose[state] += float(int(pos2) - int(pos1))
        except KeyError:
            state_counts_diagnose[state] = float(int(pos2) - int(pos1))
        
        # Write to the new BED file with translated state
        fh_out.write(f"{chrom}\t{pos1}\t{pos2}\t{state_translated}\n")
            
            
        



    
        
########################################################################################################################    

# Define the file paths
translated_file_path = "../results/translated_counts.txt"
state_file_path = "../results/state_counts.txt"

# Sort keys based on the first field
sorted_keys = sorted(state_counts_diagnose.keys(), key=lambda x: int(x[1:]))

# Open the files in write mode
with open(translated_file_path, 'w') as translated_file, open(state_file_path, 'w') as state_file:
    # Redirect output to the translated_file
    for key in sorted_keys:
        value = state_counts_diagnose[key]
        state_file.write(f"{key}\t{value/sum(state_counts_diagnose.values())}\t{state_counts_relapse[key] / sum(state_counts_relapse.values())}\n")
        
    for key, value in translated_state_counts_diagnose.items():
        translated_file.write(f"{key}\t{value/sum(translated_state_counts_diagnose.values())}\t{translated_state_counts_relapse[key] / sum(translated_state_counts_relapse.values())}\n")

        
#import matplotlib.pyplot as plt


def plot_histograms_for_key_pairs(dict1, dict2, num_rows, num_cols):
    keys = set(dict1.keys()).union(dict2.keys())
    num_keys = len(keys)

    fig, axs = plt.subplots(num_rows, num_cols, figsize=(6*num_cols, 4*num_rows))

    # Iterate over keys
    for i, key in enumerate(keys):
        row = i // num_cols
        col = i % num_cols

        value1 = dict1.get(key, 0) / sum(dict1.values())  # Get value from dict1, default to 0 if key not found
        value2 = dict2.get(key, 0) / sum(dict2.values()) # Get value from dict2, default to 0 if key not found
        
        # Plot bar chart for the current key
        axs[row, col].bar(['D', 'R'], [value1, value2], color=['blue', 'orange'])
        axs[row, col].set_ylabel(f'G. prop {key}')
        #axs[row, col].set_title(f'Values for {key}')
        axs[row, col].grid(True)

    plt.tight_layout()
    plt.show()
# Example usage:

#plot_histograms_for_key_pairs(translated_state_counts_diagnose, translated_state_counts_relapse, 4, 4)

#plot_histograms_for_key_pairs(state_counts_diagnose, state_counts_relapse, 5, 3)

#######################################################################################################################

#Get states per 100kb (Note: with this function, if a mark falls within two different 100kb windows, it is counted +1 in both)

def aggregate_bins_from_file(file_path):
    portion_counts_by_chrom = {}
    with open(file_path, 'r') as fh:
        for line in fh:
            chrom, pos1, pos2, state = line.split()
            state = translator[state]
            start = int(pos1)
            end = int(pos2)
            for portion_start in range(start // reso * reso, end, reso):
                portion_end = portion_start + reso
                if chrom not in portion_counts_by_chrom:
                    portion_counts_by_chrom[chrom] = {}
                if portion_start not in portion_counts_by_chrom[chrom]:
                    portion_counts_by_chrom[chrom][portion_start] = {}
                portion_counts = portion_counts_by_chrom[chrom][portion_start]
                if state not in portion_counts:
                    portion_counts[state] = 1
                else:
                    portion_counts[state] += 1
    return portion_counts_by_chrom


# Example usage
file_path = '../chromhmm_outputs/chromhmm_15_total_states/BALL_total_model/diagnose_15_segments.bed'
portion_counts_by_chrom = aggregate_bins_from_file(file_path)

""" # Print the result
for chrom, portion_counts in portion_counts_by_chrom.items():
    if chrom != "10":
        continue
    print(f"Chromosome {chrom}:")
    for portion_start, state_counts in sorted(portion_counts.items()):
        portion_end = portion_start + 10000
        print(f"  Portion {portion_start}-{portion_end}:")
        for state, count in state_counts.items():
            print(f"    State {state}: {count} times") """
            
###################################################################################################################

#Now we just need to assign each 100kb bin to either open, closed or polycomb based on the marks in each bin, basically do the previous but using a dictionary that takes state and turns into open, closed polycomb

#All enhancers, promoter active and poised

#For closed: Het, Sil, None, Low signal

#For Polycomb: PCdom, PCdomA
open_closed_pol_dict = {}

for chrom, portion_counts in portion_counts_by_chrom.items():
    print(f"Chromosome {chrom}:")
    for portion_start, state_counts in sorted(portion_counts.items()):
        genome_bin = f"chr{chrom}_{int(portion_start/reso)}"
        portion_end = portion_start + reso
        found_condition = False
        for state, count in state_counts.items():
            if state in ["EA", "EPr", "PoisedP", "PCdomA", "K27ac", "Sil"]:
                #print(f"  Portion {portion_start}-{portion_end} is open")
                open_closed_pol_dict[genome_bin] = "open"
                found_condition = True
                break  # Exit the inner loop
        """ elif state in ["PCdomA", "PCdom", "PolycombRepressed"]:
                #print(f"  Portion {portion_start}-{portion_end} is polycomb")
                open_closed_pol_dict[genome_bin] = "polycomb"
                found_condition = True
                break  # Exit the inner loop """

        if not found_condition:
            #print(f"  Portion {portion_start}-{portion_end} is closed")
            open_closed_pol_dict[genome_bin] = "closed"

#print(open_closed_pol_dict.items())

file_path = '../results/genomic_state_diagnose_20kb.csv'
import csv

# Write data to CSV file
with open(file_path, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    # Write header
    writer.writerow(['Bin', 'State'])
    # Write data rows
    for key, value in open_closed_pol_dict.items():
        writer.writerow([key, value])

print(f'Data has been written to {file_path}')

def calculate_percentage(dictionary):
    total_count = len(dictionary.values())
    percentages = {'open': 0, 'closed': 0}
    for value in dictionary.values():
        percentages['open'] += value.count('open')
        percentages['closed'] += value.count('closed')
        #percentages['polycomb'] += value.count('polycomb')
    for key in percentages:
        percentages[key] = (percentages[key] / total_count) * 100
    return percentages

# Calculate percentages
percentages = calculate_percentage(open_closed_pol_dict)

# Print percentages
for category, percentage in percentages.items():
    print(f'{category}: {percentage:.2f}%')
