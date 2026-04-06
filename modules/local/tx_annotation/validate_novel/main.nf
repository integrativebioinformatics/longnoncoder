process BAMBU_VALIDATE {
    tag "Validating_BAMBU"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://lfreitasl/bambu:3.8.0':
        'docker.io/lfreitasl/bambu:3.8.0' }"

    input:
    path metadata_csv
    path counts_gene
    path counts_transcript
    path cpm_transcript
    path full_length_counts_transcript
    path unique_counts_transcript

    output:
    path "BambuOutput_counts_gene_validated.txt"                 , emit: counts_gene_validated
    path "BambuOutput_counts_transcript_validated.txt"           , emit: counts_transcript_validated
    path "BambuOutput_CPM_transcript_validated.txt"              , emit: cpm_transcript_validated
    path "BambuOutput_fullLengthCounts_transcript_validated.txt" , emit: full_length_counts_transcript_validated
    path "BambuOutput_uniqueCounts_transcript_validated.txt"     , emit: unique_counts_transcript_validated
    path "versions.yml"                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # Create the validation script
    cat > validate_bambu_files.sh << 'EOF'
    # Validate against novel transcripts metadata.
    # This script reads the '_filtered.txt' files 
    # and validates BambuTx/BambuGene IDs against the metadata CSV file.

    # --- Configuration and File Paths ---

    # Define input files
    SCRIPT_DIR="\$(pwd)"
    METADATA_FILE="\${SCRIPT_DIR}/${metadata_csv}"
    TEMP_ID_FILE="\${SCRIPT_DIR}/temp_valid_ids.txt"

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
        echo "This script processes predefined '_filtered.txt' files and validates them."
        echo "Files processed:"
        for file in "\${INPUT_FILES[@]}"; do
            echo "  - \$(basename "\$file")"
        done
        echo "Output files will have '_validated' suffix."
        exit 1
    fi

    # --- Validation Functions ---

    # Function to check if metadata file exists
    check_metadata_file() {
        if [ ! -f "\$METADATA_FILE" ] || [ ! -r "\$METADATA_FILE" ]; then
            echo "Error: Metadata file '\$METADATA_FILE' not found or is not readable."
            exit 1
        fi
    }

    # Function to extract valid IDs from metadata CSV
    extract_valid_ids() {
        # Extract qry_id (column 2) and qry_gene_id (column 4) from CSV, skip header
        # Remove quotes and filter out empty values and "-"
        awk -F',' 'NR > 1 { 
            gsub(/"/, "", \$2); gsub(/"/, "", \$4); 
            if (\$2 != "" && \$2 != "-") print \$2; 
            if (\$4 != "" && \$4 != "-") print \$4 
        }' "\$METADATA_FILE" | sort -u > "\$TEMP_ID_FILE"
    }

    # --- Processing Functions ---

    # Function to process gene count file
    process_gene_file() {
        local input_file="\$1"
        local output_file="\$2"

        echo "Processing gene file: \$(basename "\$input_file")"

        # Use awk to process the file with validation
        awk -v temp_file="\$TEMP_ID_FILE" '
        BEGIN {
            # Read valid IDs into an array
            while ((getline line < temp_file) > 0) {
                valid_ids[line] = 1
            }
            close(temp_file)
        }
        NR == 1 {
            # Always print header
            print \$0
            next
        }
        {
            # If GENEID starts with "BambuGene", validate against metadata
            if (\$1 ~ /^BambuGene/) {
                if (\$1 in valid_ids) {
                    print \$0
                }
            } else {
                # For non-BambuGene IDs, just print
                print \$0
            }
        }' "\$input_file" > "\$output_file"
    }

    # Function to process transcript files
    process_transcript_file() {
        local input_file="\$1"
        local output_file="\$2"

        echo "Processing transcript file: \$(basename "\$input_file")"

        # Use awk to process the file with validation
        awk -v temp_file="\$TEMP_ID_FILE" '
        BEGIN {
            # Read valid IDs into an array
            while ((getline line < temp_file) > 0) {
                valid_ids[line] = 1
            }
            close(temp_file)
        }
        NR == 1 {
            # Always print header
            print \$0
            next
        }
        {
            # If TXNAME starts with "BambuTx", validate against metadata
            if (\$1 ~ /^BambuTx/) {
                if (\$1 in valid_ids) {
                    print \$0
                }
            } else {
                # For non-BambuTx IDs, just print
                print \$0
            }
        }' "\$input_file" > "\$output_file"
    }

    # --- Main Logic ---

    echo "Starting validation and subsetting process..."

    # Check if metadata file exists
    check_metadata_file

    # Extract valid IDs from metadata
    echo "Extracting valid IDs from metadata..."
    extract_valid_ids

    # Check if temp file was created
    if [ ! -s "\$TEMP_ID_FILE" ]; then
        echo "Error: No valid IDs were extracted from metadata. '\$TEMP_ID_FILE' is empty."
        exit 1
    fi

    # Process each input file
    for input_file in "\${INPUT_FILES[@]}"; do
        # Check if input file exists (it should be output of script 1)
        if [ ! -f "\$input_file" ] || [ ! -r "\$input_file" ]; then
            echo "Warning: Input file '\$input_file' not found or is not readable. Skipping."
            echo "         (Did you run Script 1 first?)"
            continue
        fi
        
        # Generate output filename
        base_name=\$(basename "\$input_file" _filtered.txt)
        output_file="\${SCRIPT_DIR}/\${base_name}_validated.txt"
        
        # Prevent overwriting the input file
        if [ "\$(realpath "\$input_file")" == "\$(realpath "\$output_file")" ]; then
            echo "Error: Input and output filenames would be the same for '\$input_file'. Skipping."
            continue
        fi
        
        # Determine file type and process accordingly
        if [[ "\$input_file" == *"counts_gene_filtered.txt" ]]; then
            process_gene_file "\$input_file" "\$output_file"
        else
            process_transcript_file "\$input_file" "\$output_file"
        fi
        
        # Check if output file was created successfully
        if [ -s "\$output_file" ]; then
            echo "✓ Validated subset file '\$(basename "\$output_file")' created successfully."
        else
            if [ -f "\$output_file" ]; then
                echo "⚠ Validated subset file '\$(basename "\$output_file")' created, but it might be empty or only contain the header."
            else
                echo "✗ Error: Failed to create validated subset file '\$(basename "\$output_file")'."
            fi
        fi
        
        echo ""
    done

    # Clean up temporary file
    rm -f "\$TEMP_ID_FILE"

    echo "Validation and subsetting process completed."
    echo ""
    echo "Generated files:"
    for input_file in "\${INPUT_FILES[@]}"; do
        base_name=\$(basename "\$input_file" _filtered.txt)
        output_file="\${SCRIPT_DIR}/\${base_name}_validated.txt"
        if [ -f "\$output_file" ]; then
            echo "  - \$(basename "\$output_file")"
        fi
    done    
    EOF

    # Make the script executable
    chmod +x validate_bambu_files.sh

    # Run the validation script
    ./validate_bambu_files.sh

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n1 | sed 's/GNU Awk //; s/,.*//')
        bash: \$(bash --version | head -n1 | sed 's/GNU bash, version //; s/ .*//')
    END_VERSIONS
    """

    stub:
    """
    touch BambuOutput_counts_gene_validated.txt
    touch BambuOutput_counts_transcript_validated.txt
    touch BambuOutput_CPM_transcript_validated.txt
    touch BambuOutput_fullLengthCounts_transcript_validated.txt
    touch BambuOutput_uniqueCounts_transcript_validated.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n1 | sed 's/GNU Awk //; s/,.*//')
        bash: \$(bash --version | head -n1 | sed 's/GNU bash, version //; s/ .*//')
    END_VERSIONS
    """
}