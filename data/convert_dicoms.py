import sys
import os
import pydicom
import numpy as np
from PIL import Image

def dicom_to_png(dicom_path, output_path):
    """
    Converts a DICOM file to a high-quality PNG image.

    This function reads a DICOM file, applies appropriate windowing for contrast,
    and saves the pixel data as a PNG image in the specified output path.
    If the DICOM contains multiple frames, it processes only the first one.

    Args:
        dicom_path (str): The full path to the input DICOM file.
        output_path (str): The full path where the output PNG file will be saved.
    """
    try:
        # Read the DICOM file
        ds = pydicom.dcmread(dicom_path)

        # Get pixel data as a numpy array
        pixel_data = ds.pixel_array

        # If it's a 3D volume (multiple slices), take the first slice
        if len(pixel_data.shape) > 2:
            pixel_data = pixel_data[0]

        # Normalize pixel values for display
        # First, try to use DICOM windowing information for best contrast
        if 'WindowCenter' in ds and 'WindowWidth' in ds:
            # Handle cases where window values are lists
            window_center = ds.WindowCenter[0] if isinstance(ds.WindowCenter, pydicom.multival.MultiValue) else ds.WindowCenter
            window_width = ds.WindowWidth[0] if isinstance(ds.WindowWidth, pydicom.multival.MultiValue) else ds.WindowWidth

            min_val = window_center - (window_width / 2)
            max_val = window_center + (window_width / 2)

            # Clip and scale the pixel data
            pixel_data = np.clip(pixel_data, min_val, max_val)
            pixel_data = ((pixel_data - min_val) / (max_val - min_val)) * 255.0
        else:
            # If no windowing info, use simple min-max normalization as a fallback
            pixel_data = ((pixel_data - pixel_data.min()) / (pixel_data.max() - pixel_data.min())) * 255.0

        # Convert to an 8-bit unsigned integer array
        pixel_data = pixel_data.astype(np.uint8)

        # Create a PIL Image from the numpy array
        img = Image.fromarray(pixel_data, mode='L')  # 'L' mode is for grayscale

        # Save the image as a PNG file
        img.save(output_path, 'PNG')
        print(f"✅ Successfully converted: {dicom_path} -> {output_path}")

    except Exception as e:
        print(f"❌ Could not convert {dicom_path}. Error: {e}")

def process_directory(root_path):
    """
    Recursively walks through a directory and converts all .dcm files to .png.
    """
    print(f"Scanning directory: {root_path}")
    file_count = 0
    for dirpath, _, filenames in os.walk(root_path):
        for filename in filenames:
            # Check for .dcm extension (case-insensitive)
            if filename.lower().endswith('.dcm'):
                file_count += 1
                dicom_path = os.path.join(dirpath, filename)

                # Define the output path to be next to the original file
                png_path = os.path.splitext(dicom_path)[0] + '.png'

                dicom_to_png(dicom_path, png_path)

    if file_count == 0:
        print("No DICOM (.dcm) files found in the specified directory.")
    else:
        print(f"\nProcessing complete. Found and attempted to convert {file_count} DICOM files.")


if __name__ == "__main__":
    # Ensure a directory path is provided as a command-line argument
    if len(sys.argv) != 2:
        print("Usage: python script.py <path_to_directory>")
        sys.exit(1)

    root_directory = sys.argv[1]

    # Check if the provided path is a valid directory
    if not os.path.isdir(root_directory):
        print(f"Error: The path '{root_directory}' is not a valid directory.")
        sys.exit(1)

    process_directory(root_directory)