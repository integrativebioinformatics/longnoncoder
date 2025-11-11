process SUBSET_BAMBU_GTF {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://lfreitasl/bambu:3.8.0':
        'docker.io/lfreitasl/bambu:3.8.0' }"

    input:
    tuple val(meta), path(counts_transcript) 
    tuple val(meta2), path(full_length_counts_transcript)
    tuple val(meta3), path(unique_counts_transcript)

    output:
    tuple val(meta), path("*_validated.gtf"), emit: validated_files
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # Create the validation script
    cat > subset_bambu_gtf.sh << 'EOF'
    # Script to subset GTF file based on validated transcript IDs
# Creates separate GTF files for counts_transcript, fullLengthCounts_transcript, 
# and uniqueCounts_transcript validated files.

# --- Configuration and File Paths ---

SCRIPT_DIR="\$(pwd)"
GTF_FILE="${SCRIPT_DIR}/BambuOutput_extended_annotations.gtf"

# Input validated files
declare -a VALIDATED_FILES=(
    "${SCRIPT_DIR}/BambuOutput_counts_transcript_validated.txt"
    "${SCRIPT_DIR}/BambuOutput_fullLengthCounts_transcript_validated.txt"
    "${SCRIPT_DIR}/BambuOutput_uniqueCounts_transcript_validated.txt"
)

# Corresponding output GTF files
declare -a OUTPUT_GTF_FILES=(
    "${SCRIPT_DIR}/BambuOutput_annotations_validated.gtf"
    "${SCRIPT_DIR}/BambuOutput_fullLength_validated.gtf"
    "${SCRIPT_DIR}/BambuOutput_uniquelyMapped_validated.gtf"
)

# --- Argument Handling ---

if [ "$#" -gt 0 ]; then
    if [ "$1" = "--help" ]; then
        echo "Usage: $0"
        echo "This script creates GTF subsets based on validated transcript count files."
        echo "Input files expected:"
        for file in "${VALIDATED_FILES[@]}"; do
            echo "  - $(basename "$file")"
        done
        echo "Output GTF files will be created with corresponding names."
        exit 0
    else
        echo "Error: Unknown argument(s): $*"
        echo "Use '$0 --help' for usage information."
        exit 1
    fi
fi

# --- Validation Functions ---

# Checks if the required GTF file and at least one validated transcript file exist and are readable.
# Exits with an error message if any required file is missing or unreadable.
check_files() {
    # Check if GTF file exists
    if [ ! -f "$GTF_FILE" ] || [ ! -r "$GTF_FILE" ]; then
        echo "Error: GTF file '$GTF_FILE' not found or is not readable."
        exit 1
    fi
    
    # Check if at least one validated file exists
    local found_file=false
    for file in "${VALIDATED_FILES[@]}"; do
        if [ -f "$file" ] && [ -r "$file" ]; then
            found_file=true
            break
        fi
    done
    
    if [ "$found_file" = false ]; then
        echo "Error: No validated transcript files found. Please run validate_and_subset.sh first."
        exit 1
    fi
}

# --- Processing Functions ---

# Function to extract transcript IDs from a validated file
extract_transcript_ids() {
    local input_file="$1"
    local temp_file="$2"
    
    echo "Extracting transcript IDs from $(basename "$input_file")..."
    
    # Dynamically find TXNAME column index and extract transcript IDs
    txname_col=$(head -1 "$input_file" | awk -F'\t' '{for(i=1;i<=NF;i++) if($i=="TXNAME") print i}')
    if [ -z "$txname_col" ]; then
        echo "Error: TXNAME column not found in header of '$input_file'."
        rm -f "$temp_file"
        return 1
    fi
    awk -v col="$txname_col" 'NR > 1 { print $col }' "$input_file" | sort -u > "$temp_file"
    
    local count=$(wc -l < "$temp_file")
    echo "  Found $count unique transcript IDs"
}

# Function to subset GTF file based on transcript IDs
subset_gtf() {
    local transcript_ids_file="$1"
    local output_gtf="$2"
    
    echo "Creating GTF subset: $(basename "$output_gtf")..."

    # Use embedded awk to robustly filter GTF file based on transcript_id attribute
    awk -v ids_file="$transcript_ids_file" "
    BEGIN {
        while ((getline id < ids_file) > 0) {
            ids[id] = 1
        }
        close(ids_file)
    }
    {
        # Robustly match transcript_id attribute in GTF (handles spaces, quotes, semicolons)
        match(\$0, /transcript_id[ \t]*[\"']?([^\"'; \t]+)[\"']?[ \t]*;/, arr)
        if (arr[1] && ids[arr[1]]) {
            print \$0
        }
    }
    " "$GTF_FILE" > "$output_gtf"

    # Count lines in output (excluding comments)
    local line_count=$(grep -v "^#" "$output_gtf" | wc -l)
    echo "  GTF subset created with $line_count feature lines"
}

# --- Main Logic ---

echo "Starting GTF subsetting process..."

# Check if required files exist
check_files

# Array to keep track of temp files for cleanup
declare -a TEMP_IDS_FILES=()

# Trap to clean up temp files on exit or interruption
cleanup_temp_files() {
    for temp_file in "${TEMP_IDS_FILES[@]}"; do
        [ -f "$temp_file" ] && rm -f "$temp_file"
    done
}
trap cleanup_temp_files EXIT INT TERM

# Process each validated file
processed_files=0
for i in "${!VALIDATED_FILES[@]}"; do
    input_file="${VALIDATED_FILES[$i]}"
    output_gtf="${OUTPUT_GTF_FILES[$i]}"
    
    # Check if input file exists
    if [ ! -f "$input_file" ] || [ ! -r "$input_file" ]; then
        echo "Warning: Validated file '$input_file' not found or is not readable. Skipping."
        continue
    fi
    
    # Create temporary file for transcript IDs
    temp_ids_file="$(mktemp -p "$SCRIPT_DIR" "temp_transcript_ids_$(basename "$input_file" .txt).XXXXXX")"
    TEMP_IDS_FILES+=("$temp_ids_file")
    
    # Extract transcript IDs
    extract_transcript_ids "$input_file" "$temp_ids_file"
    
    # Check if we found any transcript IDs
    if [ ! -s "$temp_ids_file" ]; then
        echo "Warning: No transcript IDs found in '$input_file'. Skipping GTF subset creation."
        rm -f "$temp_ids_file"
        continue
    fi
    
    # Subset GTF file
    subset_gtf "$temp_ids_file" "$output_gtf"

    # Check if output file was created successfully
    if [ -f "$output_gtf" ]; then
        if [ -s "$output_gtf" ]; then
            echo "✓ GTF subset '$(basename "$output_gtf")' created successfully."
        else
            echo "⚠ GTF subset '$(basename "$output_gtf")' created, but it might be empty or only contain headers."
        fi
    else
        echo "✗ Error: Failed to create GTF subset '$(basename "$output_gtf")'."
    fi
    
    # Clean up temporary file
    rm -f "$temp_ids_file"
    # Remove temp_ids_file from TEMP_IDS_FILES to avoid double deletion in trap
    TEMP_IDS_FILES=("${TEMP_IDS_FILES[@]/$temp_ids_file}")
    processed_files=$((processed_files + 1))
    echo ""
done

if [ "$processed_files" -eq 0 ]; then
    echo "Error: No validated files were processed. Exiting."
    exit 1
fi

echo "GTF subsetting process completed."
echo ""
echo "Generated GTF files:"
for output_gtf in "${OUTPUT_GTF_FILES[@]}"; do
    if [ -f "$output_gtf" ]; then
        echo "  - $(basename "$output_gtf")"
    fi
done

echo ""
echo "Summary of files created:"
echo "1. Validated transcript count files (from validate_and_subset.sh):"
for file in "${VALIDATED_FILES[@]}"; do
    if [ -f "$file" ]; then
        echo "   - $(basename "$file")"
    fi
done
echo "2. Corresponding GTF subsets:"
for output_gtf in "${OUTPUT_GTF_FILES[@]}"; do
    if [ -f "$output_gtf" ]; then
        echo "   - $(basename "$output_gtf")"
    fi
done
EOF

    # Make the script executable
    chmod +x subset_bambu_gtf.sh

    # Run the validation script
    ./subset_bambu_gtf.sh

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk --version | head -n1 | sed 's/GNU Awk //; s/,.*//')
        bash: \$(bash --version | head -n1 | sed 's/GNU bash, version //; s/ .*//')
    END_VERSIONS

        stub:
    """
    touch BambuOutput_annotations_validated.gtf
    touch BambuOutput_fullLength_validated.gtf
    touch BambuOutput_uniquelyMapped_validated.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk --version | head -n1 | sed 's/GNU Awk //; s/,.*//' || echo "5.1.0")
        bash: \$(bash --version | head -n1 | sed 's/GNU bash, version //; s/ .*//') || echo "4.4.20")
    END_VERSIONS
    """
}