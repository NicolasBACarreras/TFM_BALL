def parse_first_file(file_path):
    mapping = {}
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            chromosome = f"chr{parts[0]}"
            start = float(parts[1])
            end = float(parts[2])
            category = parts[3]
            mapping[(chromosome, start, end)] = category
    return mapping

def check_overlap(region_start, region_end, state_start, state_end):
    return region_start <= state_end and state_start <= region_end

def classify_positions(second_file_path, mapping):
    classifications = []
    with open(second_file_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            #print(parts)
            chromosome_1 = parts[1]
            chromosome_2 = parts[7]
            
            read1_start = float(parts[5])
            read1_end = float(parts[6])
            read2_start = float(parts[11])
            read2_end = float(parts[12])
            
            classification_1 = None
            classification_2 = None
            
            for key, value in mapping.items():
                chr_key, start, end = key
                #print(key, value)
                # Check first chromosome contact
                if chr_key == chromosome_1 and check_overlap(read1_start, read1_end, start, end):
                    #print("There is an overlap")
                    if classification_1:
                        if value in ["EA", "EPr", "Sil"]:
                            classification_1 = value
                    else:
                        classification_1 = value

                # Check second chromosome contact
                if chr_key == chromosome_2 and check_overlap(read2_start, read2_end, start, end):
                    if classification_2:
                        if value in ["EA", "EPr", "Sil"]:
                            classification_2 = value
                    else:
                        classification_2 = value

            classifications.append((line.strip(), classification_1, classification_2))
    return classifications

def main():
    first_file_path = '../results/relapse_15_segments_translated.bed'
    second_file_path = '../data/preprob_notch.txt'

    mapping = parse_first_file(first_file_path)
    classifications = classify_positions(second_file_path, mapping)

    # Output classifications with original line and classifications for each position
    for classification in classifications:
        original_line, classification_1, classification_2 = classification
        print(f"{original_line}\t{classification_1}\t{classification_2}")

if __name__ == "__main__":
    main()

