process SUBSET_BAMBU_GTF {
    tag "Subsetting_BAMBU_GTF"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://lfreitasl/bambu:3.8.0':
        'docker.io/lfreitasl/bambu:3.8.0' }"

    input:
    path gtf_file
    path counts_transcript
    path full_length_counts_transcript
    path unique_counts_transcript

    output:
    path "BambuOutput_annotations_validated.gtf"    , emit: annotations_validated_gtf
    path "BambuOutput_fullLength_validated.gtf"     , emit: fullLength_validated_gtf
    path "BambuOutput_uniquelyMapped_validated.gtf" , emit: uniquelyMapped_validated_gtf
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # Create the validation script
    cat > subset_bambu_gtf.sh << 'EOF'
    # Script to subset GTF file based on validated transcript IDs

    # --- Configuration and File Paths ---

    # Fail fast and propagate errors
    set -euo pipefail

    GTF_FILE="${gtf_file}"

    declare -a VALIDATED_FILES=(
        "${counts_transcript}"
        "${full_length_counts_transcript}"
        "${unique_counts_transcript}"
    )

    # Corresponding output GTF files
    declare -a OUTPUT_GTF_FILES=(
        "BambuOutput_annotations_validated.gtf"
        "BambuOutput_fullLength_validated.gtf"
        "BambuOutput_uniquelyMapped_validated.gtf"
    )

    # --- Validation Functions ---

    check_files() {
        # Check if GTF file exists
        if [ ! -f "\$GTF_FILE" ] || [ ! -r "\$GTF_FILE" ]; then
            echo "Error: GTF file '\$GTF_FILE' not found or is not readable."
            exit 1
        fi

        # Check if at least one validated file exists
        local found_file=false
        for file in "\${VALIDATED_FILES[@]}"; do
            if [ -f "\$file" ] && [ -r "\$file" ]; then
                found_file=true
                break
            fi
        done

        if [ "\$found_file" = false ]; then
            echo "Error: No validated transcript files found."
            exit 1
        fi
    }

    # Function to extract transcript IDs from a validated file
    extract_transcript_ids() {
        local input_file="\$1"
        local temp_file="\$2"

        echo "Extracting transcript IDs from \$(basename "\$input_file")..."

        # Dynamically find TXNAME column index and extract transcript IDs
        txname_col=\$(head -1 "\$input_file" | awk -F'\t' '{for(i=1;i<=NF;i++) if(\$i=="TXNAME") print i}')
        if [ -z "\$txname_col" ]; then
            echo "Error: TXNAME column not found in header of '\$input_file'."
            rm -f "\$temp_file"
            return 1
        fi

        awk -v col="\$txname_col" 'NR > 1 { print \$col }' "\$input_file" | sort -u > "\$temp_file"

        local count=\$(wc -l < "\$temp_file")
        echo "  Found \$count unique transcript IDs"
    }

    # Function to subset gtf file based on transcript IDs extracted
    subset_gtf() {
        local transcript_ids_file="\$1"
        local output_gtf="\$2"

        echo "Creating GTF subset: \$(basename "\$output_gtf")..."

        # Use the AWK script from bin folder
        # Copy the awk script to the working directory
        cp ${projectDir}/bin/subset_gtf.awk ./
        awk -v ids_file="\$transcript_ids_file" -f "./subset_gtf.awk" "\$GTF_FILE" > "\$output_gtf"

        local line_count=\$(grep -v "^#" "\$output_gtf" | wc -l)
        echo "  GTF subset created with \$line_count feature lines"
    }

    # --- Main Logic ---

    echo "Starting GTF subsetting process..."

    # Check if required files exist
    check_files

    # Process each validated file
    processed_files=0
    for i in "\${!VALIDATED_FILES[@]}"; do
        input_file="\${VALIDATED_FILES[\$i]}"
        output_gtf="\${OUTPUT_GTF_FILES[\$i]}"

        # Check if input file exists
        if [ ! -f "\$input_file" ] || [ ! -r "\$input_file" ]; then
            echo "Warning: Validated file '\$input_file' not found or is not readable. Skipping."
            continue
        fi

        # Create temporary file for transcript IDs
        temp_ids_file="\$(mktemp -p . temp_transcript_ids.XXXXXX)"

        # Extract transcript IDs
        if ! extract_transcript_ids "\$input_file" "\$temp_ids_file"; then
            echo "Error: Failed to extract transcript IDs from '\$input_file'." >&2
            rm -f "\$temp_ids_file"
            exit 1
        fi

        # Check if we found any transcript IDs
        if [ ! -s "\$temp_ids_file" ]; then
            echo "Warning: No transcript IDs found in '\$input_file'. Skipping GTF subset creation."
            rm -f "\$temp_ids_file"
            continue
        fi

        # Subset GTF file
        subset_gtf "\$temp_ids_file" "\$output_gtf"

        # Check if output file was created successfully
        if [ -f "\$output_gtf" ] && [ -s "\$output_gtf" ]; then
            echo "✓ GTF subset '\$(basename "\$output_gtf")' created successfully."
        else
            echo "✗ Error: Failed to create GTF subset '\$(basename "\$output_gtf")'."
            rm -f "\$temp_ids_file"
            exit 1
        fi

        # Clean up temporary file immediately
        rm -f "\$temp_ids_file"
        processed_files=\$((processed_files + 1))
        echo ""
    done

    if [ "\$processed_files" -eq 0 ]; then
        echo "Error: No validated files were processed. Exiting."
        exit 1
    fi

    echo "GTF subsetting process completed."
    echo ""
    echo "Generated GTF files:"
    for output_gtf in "\${OUTPUT_GTF_FILES[@]}"; do
        if [ -f "\$output_gtf" ]; then
            echo "  - \$(basename "\$output_gtf")"
        fi
    done
    EOF

    # Make the script executable
    chmod +x subset_bambu_gtf.sh

    # Run the validation script
    ./subset_bambu_gtf.sh

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n1 | sed 's/GNU Awk //; s/,.*//')
        bash: \$(bash --version | head -n1 | sed 's/GNU bash, version //; s/ .*//')
    END_VERSIONS
    """
    
    stub:
    """
    touch BambuOutput_annotations_validated.gtf
    touch BambuOutput_fullLength_validated.gtf
    touch BambuOutput_uniquelyMapped_validated.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n1 | sed 's/GNU Awk //; s/,.*//')
        bash: \$(bash --version | head -n1 | sed 's/GNU bash, version //; s/ .*//')
    END_VERSIONS
    """
}