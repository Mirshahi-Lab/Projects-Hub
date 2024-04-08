from random import randint
from pydicom import DataElement, Dataset


def extract_pulmonary(dicom_data: DataElement) -> dict:
    pulmonary = dicom_data[0x0040, 0xa730][4][0x0040,0xa730]
    pulmonary_seq_len = len(pulmonary.value)
    pulmonary_dict = {}
    for i in range(pulmonary_seq_len):
        region = pulmonary[i][0x0040, 0xa730]
        region_name = region[0][0x0040, 0xa160].value
        region_dict = {}
        max_length = len(region.value)
        for j in range(7, max_length):
            label = region[j][0x0040, 0xa043][0][0x0008, 0x0104].value
            label = label.lower().replace(' ', '_')
            if label == 'range':
                value = region[j][0x0040, 0xa168][0].CodeMeaning.lower()
                region_dict[label] = {'value': value}
            else:
                try:
                    label_unit = region[j][0x0040, 0xa300][0][0x0040, 0x08ea][0].CodeValue
                    value = float(region[j][0x0040, 0xa300][0][0x0040, 0xa30a].value)
                except KeyError:
                    label_unit = 'NA'
                    value = 'NA'
                # Label cleanup
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
    parenchyma = dicom_data[0x0040, 0xa730][4][0x0040, 0xa730]
    parenchyma_seq_len = len(parenchyma.value)
    parenchyma_dict = {}
    for i in range(parenchyma_seq_len):
        measurements = parenchyma[i][0x0040, 0xa730]
        region = measurements[0].TextValue
        unit = measurements[7][0x0040, 0xa300][0][0x0040, 0x08ea][0].CodeValue
        value = float(measurements[7][0x0040, 0xa300][0][0x0040, 0xa30a].value)
        parenchyma_dict[region] = {'unit': unit, 'value': value}
    return parenchyma_dict


def extract_cardio(dicom_data: DataElement) -> dict:
    cardio = dicom_data[0x0040,0xa730][4][0x0040, 0xa730]  # This should always be of length 2
    heart = cardio[0]
    cac = cardio[1]
    cardio_dict = {}
    for container in [heart, cac]:
        name = container[0x0040, 0xa730][-1][0x0040, 0xa043][0][0x0008, 0x0104].value.lower().replace(' ', '_')
        value = container[0x0040, 0xa730][-1][0x0040, 0xa300][0][0x0040, 0xa30a].value
        if container == heart:
            cardio_dict[name] = {'unit': 'mL', 'value': value}
        elif container == cac:
            cardio_dict[name] = {'unit': 'mm3', 'value': value}
    return cardio_dict


def extract_aorta(dicom_data: DataElement) -> dict:
    aorta = dicom_data[0x0040, 0xa730][4][0x0040, 0xa730]  # This will always have 9 elements
    aorta_dict = {}
    for i in range(9):
        aortic_location = aorta[i][0x0040, 0xa730][-4][0x0040, 0xa168][0][0x0008, 0x0104].value
        aortic_measurement = float(aorta[i][0x0040,0xa730][-2].MeasuredValueSequence[0][0x0040,0xa30a].value)
        aorta_dict[aortic_location] = {'unit': 'mm', 'value': aortic_measurement}
    return aorta_dict


def extract_spine(dicom_data: DataElement) -> dict:
    spine = dicom_data[0x0040, 0xa730][4][0x0040, 0xa730]
    spine_len = range(len(spine.value))
    height_index_dict = {
        5: 'anterior',
        6: 'medial',
        7: 'posterior'
    }
    spine_dict = {}
    for i in spine_len:
        vert = spine[i][0x0040, 0xa730]
        vert_name = vert[0].TextValue
        vert_dict = {}
        if vert_name != 'Spine Applied Range':
            for idx, label in height_index_dict.items():
                height = vert[idx][0x0040, 0xa300][0][0x0040, 0xa30a].value
                med_range = vert[idx][0x0040, 0xa730][1][0x0040, 0xa168][0].CodeMeaning.lower()
                vert_dict[f'{label}_height'] = {'unit': 'mm', 'value': height, 'range': med_range}
            avg_hu = vert[8][0x0040, 0xa300][0][0x0040, 0xa30a].value
            vert_dict['attenuation'] = {'unit': 'HU', 'value': avg_hu}
            spine_dict[vert_name] = vert_dict
        else:
            spine_range = vert[4][0x0040, 0xa168][0].CodeMeaning.lower()
            spine_dict['spine_med_range'] = spine_range
    return spine_dict

def extract_lesions(dicom_data: DataElement) -> dict:
    lesions = dicom_data[0x0040, 0xa730][4][0x0040, 0xa730]
    num_lesions = len(lesions.value)
    lesions_dict = {}
    for i in range(num_lesions):
        lesion_data = lesions[i][0x0040, 0xa730]
        lesion_name = lesion_data[0].TextValue
        lesion_location = lesion_data[3][0x0040, 0xa730][0][0x0040, 0xa168][0].CodeMeaning.lower()
        long_axis_mm = float(lesion_data[6][0x0040, 0xa300][0][0x0040, 0xa30a].value)
        mesh_3d_diam = float(lesion_data[7][0x0040, 0xa300][0][0x0040, 0xa30a].value)
        short_axis_mm = float(lesion_data[8][0x0040, 0xa300][0][0x0040, 0xa30a].value)
        avg_2d_diameter = float(lesion_data[9][0x0040, 0xa300][0][0x0040, 0xa30a].value)
        lesion_volume = float(lesion_data[10][0x0040, 0xa300][0][0x0040, 0xa30a].value)
        lesions_dict[lesion_name] = {
            'location': lesion_location,
            'diameters': {
                'unit': 'mm',
                'long_axis': long_axis_mm,
                'short_axis': short_axis_mm,
                'avg_2d': avg_2d_diameter,
                '3d_mesh': mesh_3d_diam
            },
            'volume': lesion_volume
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
        'AcquisitionDate': {
            'default': 'YYYY-MM-DD',
            'tag_name': 'scan_date'
        },
    }
    identfiers = {}
    for tag, defaults in identifier_dict.items():
        output = getattr(dicom, tag, defaults['default'])
        if output == '':
            output = defaults['default']
        identfiers[defaults['tag_name']] = output
    return identfiers
