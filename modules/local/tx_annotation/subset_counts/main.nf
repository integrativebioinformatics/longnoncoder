process SUBSET_BAMBU_COUNTS {
    tag "Subsetting_Counts"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://lfreitasl/bambu:3.8.0':
        'docker.io/lfreitasl/bambu:3.8.0' }"

    input:
    path counts_gene
    path counts_transcript
    path cpm_transcript
    path full_length_counts_transcript
    path unique_counts_transcript

    output:
    path "BambuOutput_counts_gene_filtered.txt"                 , emit: counts_gene_filtered
    path "BambuOutput_counts_transcript_filtered.txt"           , emit: counts_transcript_filtered
    path "BambuOutput_CPM_transcript_filtered.txt"              , emit: cpm_transcript_filtered
    path "BambuOutput_fullLengthCounts_transcript_filtered.txt" , emit: full_length_counts_transcript_filtered
    path "BambuOutput_uniqueCounts_transcript_filtered.txt"     , emit: unique_counts_transcript_filtered
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    
    # Create the validation script
    cat > subset_bambu_counts.sh << 'EOF'
    # Filter out zero-count rows.
    # This script reads the original count/CPM files and removes any
    # gene or transcript where all sample columns are 0.

    # --- Configuration and File Paths ---

    # Define input files
    SCRIPT_DIR="\$(pwd)"

    # Input files to process
    declare -a INPUT_FILES=(
        "\${SCRIPT_DIR}/${counts_gene}"
        "\${SCRIPT_DIR}/${counts_transcript}"
        "\${SCRIPT_DIR}/${cpm_transcript}"
        "\${SCRIPT_DIR}/${full_length_counts_transcript}"
        "\${SCRIPT_DIR}/${unique_counts_transcript}"
    )

    # --- Argument Handling ---

    if [ "\$#" -ne 0 ]; then
        echo "Usage: \$0"
        echo "This script filters predefined input files for zero-count rows."
        echo "Files processed:"
        for file in "\${INPUT_FILES[@]}"; do
            echo "  - \$(basename "\$file")"
        done
        echo "Output files will have '_filtered' suffix."
        exit 1
    fi

    # --- Processing Functions ---

    # Function to process gene count file (checks cols 2-NF)
    process_gene_file() {
        local input_file="\$1"
        local output_file="\$2"

        echo "Processing gene file: \$(basename "\$input_file")"

        # Use awk to filter rows with all-zero counts
        awk '
        NR == 1 {
            # Always print header
            print \$0
            next
        }
        {
            # Check if any of the count columns are not zero
            has_counts = 0
            for (i = 2; i <= NF; i++) {
                if (\$i != 0) {
                    has_counts = 1
                    break
                }
            }

            if (has_counts) {
                print \$0
            }
        }' "\$input_file" > "\$output_file"
    }

    # Function to process transcript files 
    process_transcript_file() {
        local input_file="\$1"
        local output_file="\$2"

        echo "Processing transcript file: \$(basename "\$input_file")"

        # Use awk to filter rows with all-zero counts
        awk '
        NR == 1 {
            # Always print header
            print \$0
            next
        }
        {
            # Check if any of the count columns are not zero
            has_counts = 0
            for (i = 3; i <= NF; i++) {
                if (\$i != 0) {
                    has_counts = 1
                    break
                }
            }

            if (has_counts) {
                print \$0
            }
        }' "\$input_file" > "\$output_file"
    }

    # --- Main Logic ---

    echo "Starting zero-count filtering process..."

    # Process each input file
    for input_file in "\${INPUT_FILES[@]}"; do
        # Check if input file exists
        if [ ! -f "\$input_file" ] || [ ! -r "\$input_file" ]; then
            echo "Warning: Input file '\$input_file' not found or is not readable. Skipping."
            continue
        fi

        # Generate output filename
        base_name=\$(basename "\$input_file" .txt)
        output_file="\${SCRIPT_DIR}/\${base_name}_filtered.txt"

        # Prevent overwriting the input file
        if [ "\$(realpath "\$input_file")" == "\$(realpath "\$output_file")" ]; then
            echo "Error: Input and output filenames would be the same for '\$input_file'. Skipping."
            continue
        fi

        # Determine file type and process accordingly
        if [[ "\$input_file" == *"counts_gene.txt" ]]; then
            process_gene_file "\$input_file" "\$output_file"
        else
            process_transcript_file "\$input_file" "\$output_file"
        fi

        # Check if output file was created successfully
        if [ -s "\$output_file" ]; then
            echo "✓ Filtered file '\$(basename "\$output_file")' created successfully."
        else
            if [ -f "\$output_file" ]; then
                echo "⚠ Filtered file '\$(basename "\$output_file")' created, but it might be empty or only contain the header."
            else
                echo "✗ Error: Failed to create filtered file '\$(basename "\$output_file")'."
            fi
        fi

        echo ""
    done

    echo "Zero-count filtering process completed."
    echo ""
    echo "Generated files:"
    for input_file in "\${INPUT_FILES[@]}"; do
        base_name=\$(basename "\$input_file" .txt)
        output_file="\${SCRIPT_DIR}/\${base_name}_filtered.txt"
        if [ -f "\$output_file" ]; then
            echo "  - \$(basename "\$output_file")"
        fi
    done
    EOF

    # Make the script executable
    chmod +x subset_bambu_counts.sh

    # Run the validation script
    ./subset_bambu_counts.sh

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n1 | sed 's/GNU Awk //; s/,.*//')
        bash: \$(bash --version | head -n1 | sed 's/GNU bash, version //; s/ .*//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    touch BambuOutput_counts_gene_filtered.txt
    touch BambuOutput_counts_transcript_filtered.txt
    touch BambuOutput_CPM_transcript_filtered.txt
    touch BambuOutput_fullLengthCounts_transcript_filtered.txt
    touch BambuOutput_uniqueCounts_transcript_filtered.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n1 | sed 's/GNU Awk //; s/,.*//')
        bash: \$(bash --version | head -n1 | sed 's/GNU bash, version //; s/ .*//')
    END_VERSIONS
    """
}