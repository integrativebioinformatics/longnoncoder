#!/bin/bash

# Quick test using stub mode (no real computation)
echo "🧪 Running quick test with stub mode..."

nextflow run test_workflow.nf -profile conda -stub-run

if [ $? -eq 0 ]; then
    echo "✅ Stub test passed! The module structure is correct."
    echo "💡 Now you can run the full test with: ./run_test.sh"
else
    echo "❌ Stub test failed! Check your module configuration."
fi