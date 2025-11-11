import sys
import os
import pydicom

# A list of common DICOM tags that contain Protected Health Information (PHI).
# This list is not exhaustive and may need to be expanded depending on the
# specific DICOM files you are working with.
TAGS_TO_ANONYMIZE = {
    # === Patient Information ===
    'PatientName':                (0x0010, 0x0010),
    'PatientID':                  (0x0010, 0x0020),
    'PatientBirthDate':           (0x0010, 0x0030),
    'PatientSex':                 (0x0010, 0x0040),
    'PatientAge':                 (0x0010, 0x1010),
    'PatientAddress':             (0x0010, 0x1040),
    'PatientTelephoneNumbers':    (0x0010, 0x2154),
    'OtherPatientIDs':            (0x0010, 0x1000),
    'OtherPatientNames':          (0x0010, 0x1001),
    'EthnicGroup':                (0x0010, 0x2160),
    'PatientComments':            (0x0010, 0x4000),

    # === Study/Visit Information ===
    'InstitutionName':            (0x0008, 0x0080),
    'InstitutionAddress':         (0x0008, 0x0081),
    'ReferringPhysicianName':     (0x0008, 0x0090),
    'PhysiciansOfRecord':         (0x0008, 0x1048),
    'PerformingPhysicianName':    (0x0008, 0x1050),
    'NameOfPhysiciansReadingStudy': (0x0008, 0x1060),
    'OperatorsName':              (0x0008, 0x1070),
    'AdmittingDiagnosesDescription': (0x0008, 0x1080),
    'AdditionalPatientHistory':   (0x0010, 0x21B0),
    'StudyID':                    (0x0020, 0x0010),
    'AccessionNumber':            (0x0008, 0x0050),
}


def anonymize_dicom_file(filepath):
    """
    Anonymizes a single DICOM file by clearing PHI in specified tags.
    This function overwrites the original file.

    Args:
        filepath (str): The full path to the DICOM file.
    """
    try:
        # Read the DICOM file
        dataset = pydicom.dcmread(filepath)

        # A flag to check if any modifications were made
        modified = False

        # Iterate through the dictionary of tags to anonymize
        for name, tag in TAGS_TO_ANONYMIZE.items():
            if tag in dataset:
                modified = True
                # Overwrite the tag value based on its Value Representation (VR)
                if dataset[tag].VR == "PN":  # Person Name
                    dataset[tag].value = "ANONYMIZED"
                elif dataset[tag].VR in ["DA", "DT", "TM"]:  # Date, DateTime, Time
                    dataset[tag].value = ""
                else:
                    dataset[tag].value = "ANONYMIZED"

        # If any changes were made, overwrite the original file
        if modified:
            dataset.save_as(filepath)
            print(f"✅ Anonymized and saved: {filepath}")
        else:
            print(f"ℹ️ No PHI tags found to anonymize in: {filepath}")

    except Exception as e:
        print(f"❌ Could not process {filepath}. Error: {e}")


def process_directory(root_path):
    """
    Recursively walks through a directory and anonymizes all .dcm files.
    """
    print(f"\nScanning directory: {root_path}")
    file_count = 0
    for dirpath, _, filenames in os.walk(root_path):
        for filename in filenames:
            if filename.lower().endswith('.dcm'):
                file_count += 1
                full_path = os.path.join(dirpath, filename)
                anonymize_dicom_file(full_path)

    if file_count == 0:
        print("No DICOM (.dcm) files found in the specified directory.")
    else:
        print(f"\nProcessing complete. Found and attempted to anonymize {file_count} DICOM files.")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python anonymize_dicoms.py <path_to_directory>")
        sys.exit(1)

    root_directory = sys.argv[1]

    if not os.path.isdir(root_directory):
        print(f"Error: The path '{root_directory}' is not a valid directory.")
        sys.exit(1)

    # --- Safety Confirmation ---
    print("====================================================================")
    print("                         !!! WARNING !!!")
    print("====================================================================")
    print("This script will PERMANENTLY MODIFY and OVERWRITE DICOM files.")
    print("It is designed to remove Personal Health Information (PHI) from the")
    print("file headers to anonymize the data.")
    print("\n>>> PLEASE MAKE SURE YOU HAVE A BACKUP OF YOUR DATA BEFORE PROCEEDING. <<<")
    print("====================================================================")

    # Require the user to explicitly type 'yes' to proceed
    response = input(f"Are you sure you want to anonymize all .dcm files in '{root_directory}'? (Type 'yes' to confirm): ")

    if response.strip().lower() == 'yes':
        process_directory(root_directory)
    else:
        print("\nOperation cancelled by user. No files were changed.")
        sys.exit(0)