#!/bin/bash
# You can use it as ./copy_structure.sh --source_path /path/to/source --destination_path /path/to/target

# Function to display usage information
usage() {
    echo "Usage: $0 --source_path <source-dir> --destination_path <target-dir>"
    exit 1
}

# Initialize variables
SOURCE_DIR=""
TARGET_DIR=""

# Parse command-line options
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --source_path) SOURCE_DIR="$2"; shift ;;
        --destination_path) TARGET_DIR="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; usage ;;
    esac
    shift
done

# Check if both source and target directories have been set
if [ -z "${SOURCE_DIR}" ] || [ -z "${TARGET_DIR}" ]; then
    echo "Both source and destination paths must be provided."
    usage
fi

# Create the target directory if it doesn't already exist
mkdir -p "${TARGET_DIR}"

# Function to copy specific files if they exist
copy_files() {
    local src_dir=$1
    local dest_dir=$2
    for file in POSCAR INCAR OUTCAR OSZICAR POTCAR KPOINTS; do
        # Check if file exists and copy
        if [ -f "${src_dir}/${file}" ]; then
            cp "${src_dir}/${file}" "${dest_dir}"
        fi
    done
}

# Export the function so it can be used in subshells
export -f copy_files

# Find all directories in the source directory, create corresponding
# directories in the target directory, and copy specific files
find "${SOURCE_DIR}" -type d -print0 | while IFS= read -r -d $'\0' dir; do
    # Create a corresponding directory in the target location
    target_dir="${TARGET_DIR}/${dir#${SOURCE_DIR}/}"
    mkdir -p "${target_dir}"
    # Call function to copy files
    copy_files "${dir}" "${target_dir}"
done
