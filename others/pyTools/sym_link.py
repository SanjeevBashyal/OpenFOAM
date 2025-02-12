import os
import os.path
import shutil

def create_symlinks(source_folder, target_folder, subfolders_to_include):
    # Convert paths to absolute paths for accurate comparison
    source_folder = os.path.abspath(source_folder)
    target_folder = os.path.abspath(target_folder)

    # Clear the target folder if it exists
    if os.path.exists(target_folder):
        print(f"Clearing target folder: {target_folder}")
        shutil.rmtree(target_folder)  # Delete the folder and its contents

    # Recreate the target folder
    os.makedirs(target_folder)

    # Convert subfolders_to_include to absolute paths
    subfolders_to_include = [os.path.abspath(os.path.join(source_folder, subfolder)) for subfolder in subfolders_to_include]

    # Walk through the source folder
    for root, dirs, files in os.walk(source_folder):
        # Skip the target folder
        if os.path.abspath(root) == target_folder:
            continue

        # Check if the current directory is a child of any directory in subfolders_to_include
        is_child = any(os.path.commonpath([root, subfolder]) == subfolder for subfolder in subfolders_to_include)
        if not is_child:
            continue  # Skip this directory

        for file in files:
            # Get the full path of the source file
            source_file_path = os.path.join(root, file)
            
            # Use only the file name for the symlink
            symlink_name = file
            symlink_path = os.path.join(target_folder, symlink_name)
            
            # Handle conflicts if a file with the same name already exists
            if os.path.exists(symlink_path):
                print(f"Conflict: Symlink already exists for {symlink_name}. Skipping {source_file_path}")
                continue
            
            # Create the symbolic link
            os.symlink(source_file_path, symlink_path)
            print(f"Created symlink: {symlink_path} -> {source_file_path}")

# Example usage
source_folder = '/usr/lib/openfoam/openfoam2312/bashyal/Utilities/'
target_folder = '/usr/lib/openfoam/openfoam2312/bashyal/Utilities/lnInclude/'
subfolders_to_include = ['includeHeaders']

create_symlinks(source_folder, target_folder, subfolders_to_include)