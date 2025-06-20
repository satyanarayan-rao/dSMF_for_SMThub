max_windows = 'data_from_ReMap/sample_max_peaks_windows.fa'
windows = open(max_windows)
header = ""

for line in windows:
    line = line.strip()  # Ensure line is stripped of whitespace

    if line.startswith(">"):
        header = line  # Set header to the current line
    elif "CACGTG" in line or "cacgtg" in line:  # Check for motif
        header += "_wmotif\n" + line + "\n"
        print(header)  # Print after processing this sequence with motif
    else:
        header += "_womotif\n" + line + "\n"
        print(header)  # Print after processing this sequence without motif

windows.close()  # Close the file when done
