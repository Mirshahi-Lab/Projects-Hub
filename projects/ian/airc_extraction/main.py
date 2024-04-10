import argparse
import os
import json
from glob import glob
from pathlib import Path
from pydicom import dcmread
from pydicom.errors import InvalidDicomError

from extraction import (
    extract_identifiers,
    extract_pulmonary,
    extract_parenchyma,
    extract_cardio,
    extract_aorta,
    extract_spine,
    extract_lesions
)


# Lots of magic tags and numbers from the dicoms - it is what it is
def main(args):
    # Check that the input directory is real
    if not os.path.exists(args.input_folder):
        raise FileNotFoundError(f'{args.input_folder} could not be found. Please check your input')
    input_dicoms = glob(f'{args.input_folder}/*')
    # Check that there are actually dicoms to extract
    if not input_dicoms:
        raise ValueError(f'{args.input_folder} appears to be empty, please check that the AIRC output is there.')
    if args.output_file is not None:
        base_folder = os.path.dirname(args.output_file)
        if not os.path.exists(base_folder):
            os.makedirs(base_folder)
    else:
        args.output_file = f'{args.input_folder}/summary.json'
    # Loop through the dicoms
    structured_report = []
    for dicom in input_dicoms:
        try:
            dcm = dcmread(dicom)
            structured_report.append(dcm)
        except InvalidDicomError:
            print(f'{dicom} is not a valid dicom file, will be skipped in structured report')
    print(f'Loaded {len(structured_report)} dicom files from the input folder.')
    # Take the first instance from the structured report and load the identifying header data as the start of the output
    output_dict = extract_identifiers(structured_report[0])
    # loop through the structured report, matching cases to appropriate
    for report in structured_report:
        # Will need to implement error handling logic here in the future once we get non-dummy data
        # This extracts the name of the dicom report type
        # 
        # report_name = report[0x0040, 0xa730][0][0x0040, 0xa043][0][0x0008, 0x0104].value # this is from the example
        report_name = report[0x0040, 0xa043][0].CodeMeaning
        match report_name:
            case 'AI-Rad CT Pulmonary Density':
                pulmonary_data = extract_pulmonary(report)
                output_dict['pulmonary'] = pulmonary_data
            case 'AI-Rad CT Lung Parenchyma':
                parencyhma_data = extract_parenchyma(report)
                output_dict['parenchyma'] = parencyhma_data
            case 'AI-Rad CT Cardio':
                cardio_data = extract_cardio(report)
                output_dict['cardio'] = cardio_data
            case 'AI-Rad CT Vascular Aorta':
                aorta_data = extract_aorta(report)
                output_dict['aorta'] = aorta_data
            case 'AI-Rad CT Spine':
                spine_data = extract_spine(report)
                output_dict['spine'] = spine_data
            case 'AI-Rad CT Lung Lesion':
                lesion_data = extract_lesions(report)
                output_dict['lesions'] = lesion_data
    # Write the output_dict to the output json file:
    with open(args.output_file, 'w') as f:
        json.dump(output_dict, f, indent=4)
    print(f'Summary output to {args.output_file}!')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract data from Siemens AI-Rad Companion Structured reports and '
                                                 'save them to json format')
    parser.add_argument('-i', '--input-folder', type=Path, required=True, help='Path to folder containing the structured '
                                                                        'report as dicom files')
    parser.add_argument('-o', '--output-file', type=Path, default=None, help='Path to output file, defaults to '
                                                                             'input_folder/summary.json')
    args = parser.parse_args()
    main(args)
