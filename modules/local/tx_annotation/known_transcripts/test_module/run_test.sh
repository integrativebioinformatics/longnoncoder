#!/bin/bash

# Test script for transcript annotation module

echo "🧪 Starting transcript annotation module test..."
echo "================================================"

# Clean up previous runs
if [ -d "results" ]; then
    echo "🧹 Cleaning up previous results..."
    rm -rf results
fi

if [ -d "work" ]; then
    echo "🧹 Cleaning up work directory..."
    rm -rf work
fi

# Run the test with conda profile
echo "🚀 Running test with conda profile..."
nextflow run test_workflow.nf -profile singularity

# Check if the run was successful
if [ $? -eq 0 ]; then
    echo "✅ Test completed successfully!"
    echo ""
    echo "📊 Results summary:"
    echo "==================="
    
    if [ -d "results/transcript_annotation" ]; then
        echo "📁 Output files generated:"
        ls -la results/transcript_annotation/
        echo ""
        
        # Check key output files
        key_files=(
            "annotated_transcriptome_metadata.csv"
            "annotated_lncRNAs_metadata.csv" 
            "annotated_protein-coding_metadata.csv"
            "bambu_annotated_transcriptome.gtf"
        )
        
        echo "🔍 Checking key output files:"
        for file in "${key_files[@]}"; do
            if [ -f "results/transcript_annotation/$file" ]; then
                size=$(wc -l < "results/transcript_annotation/$file")
                echo "  ✅ $file ($size lines)"
            else
                echo "  ❌ $file (missing)"
            fi
        done
    else
        echo "❌ Results directory not found!"
    fi
    
else
    echo "❌ Test failed!"
    echo "Check the .nextflow.log file for details"
    exit 1
fi

echo ""
echo "🎉 Test completed! Check the results directory for outputs."
