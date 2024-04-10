from datetime import datetime
from random import randint
from pydicom import DataElement, Dataset


def extract_pulmonary(dicom_data: DataElement) -> dict:
    pulmonary = dicom_data[0x0040, 0xa730][-1][0x0040,0xa730]
    pulmonary_seq_len = len(pulmonary.value)
    pulmonary_dict = {}
    for i in range(pulmonary_seq_len):
        region_name = pulmonary[i][0x0040, 0xa043][0].CodeMeaning.lower().replace(' ', '_')
        region_measurements = pulmonary[i][0x0040, 0xa730]
        num_measurements = len(region_measurements.value)
        region_dict = {}
        for j in range(num_measurements):
            label = region_measurements[j][0x0040, 0xa043][0].CodeMeaning
            label = label.lower().replace(' ', '_')
            if label == 'range':
                value = region[j][0x0040, 0xa168][0].CodeMeaning.lower()
                region_dict[label] = {'value': value}
            else:
                try:
                    label_unit = region_measurements[j][0x0040, 0xa300][0][0x0040, 0x08ea][0].CodeValue
                    value = float(region_measurements[j][0x0040, 0xa300][0][0x0040, 0xa30a].value)

                except KeyError:
                    label_unit = 'NA'
                    value = 'NA'
                if label_unit == '{Number}':
                    label_unit = 'unitless'
                elif label_unit == "[hnsf'U]":
                    label_unit = 'HU'
                elif label_unit == '%{vol}':
                    label_unit = '%{total_volume}'
                region_dict[label] = {'unit': label_unit, 'value': value}
        pulmonary_dict[region_name] = region_dict
    return pulmonary_dict


def extract_parenchyma(dicom_data: DataElement) -> dict:
    parenchyma = dicom_data[0x0040, 0xa730][-1][0x0040, 0xa730]
    parenchyma_seq_len = len(parenchyma.value)
    parenchyma_dict = {}
    for i in range(parenchyma_seq_len):
        measurements = parenchyma[i]
        region = measurements[0x0040, 0xa043][0].CodeMeaning.lower().replace(' ', '_')
        unit = measurements[0x0040, 0xa730][0][0x0040, 0xa300][0][0x0040, 0x08ea][0].CodeValue
        value = float(measurements[0x0040, 0xa730][0][0x0040, 0xa300][0][0x0040, 0xa30a].value)
        parenchyma_dict[region] = {'unit': unit, 'value': value}
    return parenchyma_dict


def extract_cardio(dicom_data: DataElement) -> dict:
    cardio = dicom_data[0x0040,0xa730][-1][0x0040, 0xa730]  # This should always be of length 2
    cardio_len = len(cardio.value)
    cardio_dict = {}
    # for container in [heart, cac]:
    #     name = container[0x0040, 0xa730][-1][0x0040, 0xa043][0][0x0008, 0x0104].value.lower().replace(' ', '_')
    #     value = container[0x0040, 0xa730][-1][0x0040, 0xa300][0][0x0040, 0xa30a].value
    #     if container == heart:
    #         cardio_dict[name] = {'unit': 'mL', 'value': value}
    #     elif container == cac:
    #         cardio_dict[name] = {'unit': 'mm3', 'value': value}
    for idx in range(cardio_len):
        name = cardio[idx][0x0040, 0xa043][0].CodeMeaning.lower().replace(' ', '_')
        value = cardio[idx][0x0040, 0xa730][-1][0x0040, 0xa300][0][0x0040, 0xa30a].value
        if name == 'heart':
            cardio_dict[name] = {'unit': 'mL', 'value': value}
        elif name == 'total_coronary_calcium':
            cardio_dict[name] = {'unit': 'mm3', 'value': value}
    return cardio_dict


def extract_aorta(dicom_data: DataElement) -> dict:
    aorta = dicom_data[0x0040, 0xa730][-1][0x0040, 0xa730]  # This will always have 9 elements
    aorta_dict = {}
    for i in range(9):
        location = aorta[i][0x0040, 0xa043][0].CodeMeaning.lower().replace(' ', '_')
        diameter = float(aorta[i][0x0040, 0xa730][0][0x0040, 0xa300][0][0x0040, 0xa30a].value)
        aorta_dict[location] = {'unit': 'mm', 'value': diameter}
    return aorta_dict


def extract_spine(dicom_data: DataElement) -> dict:
    spine = dicom_data[0x0040, 0xa730][-1][0x0040, 0xa730]
    spine_len = range(len(spine.value))
    height_index_dict = {
        0: 'anterior',
        1: 'medial',
        2: 'posterior'
    }
    spine_dict = {}
    spine_name_dict = {
        'First thoracic vertebra': 'T1',
        'Second thoracic vertebra': 'T2',
        'Third thoracic vertebra': 'T3',
        'Fourth thoracic vertebra': 'T4',
        'Fifth thoracic vertebra': 'T5',
        'Sixth thoracic vertebra': 'T6',
        'Seventh thoracic vertebra': 'T7',
        'Eighth thoracic vertebra': 'T8',
        'Ninth thoracic vertebra': 'T9',
        'Tenth thoracic vertebra': 'T10',
        'Eleventh thoracic vertebra': 'T11',
        'Twelfth thoracic vertebra': 'T12'
    }
    for i in spine_len:
        vert_dict = {}
        name = spine[i][0x0040, 0xa043][0].CodeMeaning
        if name != 'Spine Applied Range':
            vert = spine[i][0x0040, 0xa730]
            vert_name = spine_name_dict.get(name)
            for idx, label in height_index_dict.items():
                height = vert[idx][0x0040, 0xa730][0][0x0040, 0xa300][0][0x0040, 0xa30a].value
                med_range = vert[idx][0x0040, 0xa730][1][0x0040, 0xa168][0].CodeMeaning.lower()
                vert_dict[f'{label}_height'] = {'unit': 'mm', 'value': height, 'range': med_range}
            avg_hu = vert[-1][0x0040, 0xa300][0][0x0040, 0xa30a].value
            vert_dict['attenuation'] = {'unit': 'HU', 'value': avg_hu}
            spine_dict[vert_name] = vert_dict
        else:
            vert = spine[i]
            spine_range = vert[0x0040, 0xa168][0].CodeMeaning.lower()
            spine_dict['spine_med_range'] = spine_range
    return spine_dict


def extract_lesions(dicom_data: DataElement) -> dict:
    lesions = dicom_data[0x0040, 0xa730][-1][0x0040, 0xa730]
    num_lesions = len(lesions.value)
    lesions_dict = {}
    for i in range(num_lesions):
        if lesions[i][0x0040, 0xa043][0].CodeMeaning != 'Tumor Burden':
            lesion_data = lesions[i][0x0040, 0xa730]
            lesion_name = lesion_data[0].TextValue
            lesion_location = lesion_data[4][0x0040, 0xa168][0].CodeMeaning.lower()
            recist_mm = float(lesion_data[1][0x0040, 0xa300][0][0x0040, 0xa30a].value)
            mesh_3d_diam = float(lesion_data[2][0x0040, 0xa300][0][0x0040, 0xa30a].value)
            lesion_volume = float(lesion_data[3][0x0040, 0xa300][0][0x0040, 0xa30a].value)
            review_status = lesion_data[5][0x0040, 0xa160].value.lower()
            if 'unreviewed' in review_status:
                review_status = 'unreviewed'
            lesions_dict[lesion_name] = {
                'location': lesion_location,
                'diameters': {
                    'unit': 'mm',
                    'recist': recist_mm,
                    '3d_mesh': mesh_3d_diam
                },
                'volume': {'unit': 'mm3', 'value': lesion_volume},
                'review_status': review_status
            }
    return lesions_dict


def extract_identifiers(dicom: Dataset) -> dict:
    identifier_dict = {
        'PatientID': {
            'default': f'FakeMRN{randint(0, 1e8)}',
            'tag_name': 'mrn'
        },
        'AccessionNumber': {
            'default': f'FakeAcc{randint(0, 1e8)}',
            'tag_name': 'accession'
        },
        'SeriesDescription': {
            'default': 'NoSeriesFound',
            'tag_name': 'series_name'
        },
        'SequenceName': {
            'default': 'NoSequenceFound',
            'tag_name': 'sequence_name'
        },
        'ProtocolName': {
            'default': 'NoProtocolFound',
            'tag_name': 'protocol_name'
        },
        'StudyDate': {
            'default': 'YYYY-MM-DD',
            'tag_name': 'scan_date'
        },
    }
    identfiers = {}
    for tag, defaults in identifier_dict.items():
        output = getattr(dicom, tag, defaults['default'])
        if output == '':
            output = defaults['default']
        elif tag == 'StudyDate':
            output = datetime.strptime(output, '%Y%m%d').date()
            output = output.strftime('%Y-%m-%d')
        identfiers[defaults['tag_name']] = output
    return identfiers
