from inquant_tools import *

path_to_input_files = 'Datafiles/mixed_proteome/'

predictions_file = path_to_input_files + 'mixed.csv'
mzml_files = [
    path_to_input_files + 'mzml/20240221_MK_MEKONG_KK_FAIMS_1CV_Evo_40SPD_Whisper_1475_mixed_proteome_A1.mzML',
    path_to_input_files + 'mzml/20240221_MK_MEKONG_KK_FAIMS_1CV_Evo_40SPD_Whisper_1475_mixed_proteome_A2.mzML',
    path_to_input_files + 'mzml/20240221_MK_MEKONG_KK_FAIMS_1CV_Evo_40SPD_Whisper_1475_mixed_proteome_A3.mzML',
    path_to_input_files + 'mzml/20240221_MK_MEKONG_KK_FAIMS_1CV_Evo_40SPD_Whisper_1475_mixed_proteome_B1.mzML',
    path_to_input_files + 'mzml/20240221_MK_MEKONG_KK_FAIMS_1CV_Evo_40SPD_Whisper_1475_mixed_proteome_B2.mzML',
    path_to_input_files + 'mzml/20240221_MK_MEKONG_KK_FAIMS_1CV_Evo_40SPD_Whisper_1475_mixed_proteome_B3.mzML'
]

proteome_file = path_to_input_files + 'mixed_human_yeast_ecoli.fasta'

output_folder = 'Results/mixed_proteome/'

# Getting the tolerance from the calibration
with INQuant(predictions_file=predictions_file,
             mzml_file_list=mzml_files,
             script_purpose='Mixed proteome calibration',
             experiment_name='mixed_proteome_calibration',
             confidence_filter=0.0,
             output_file_path=output_folder) as INQ:
    
    # Load predictions using INQuant
    predictions_dict = INQ.load_predictions()

    calibration = Calibration(predictions_dict)

    for mzml_id, mzml_values in INQ.mzml_file_dict.items():
        path = mzml_values['file path']
        name = mzml_values['file name']
        calibration.load_mzml(path)
        calibration.calibration(mzml_id, set_training_seed=42)
        mzml_max_tolerance = calibration.ppm_tolerances[mzml_id]
        INQ.variable_status.append(f"- Maximum ppm tolerance for {mzml_id}: {mzml_max_tolerance}\n")

    new_tolerance = calibration.get_ppm_tolerance()
    INQ.variable_status.append(f"- Maximum ppm tolerance across all experiments: {new_tolerance}\n")

# Running with MBR
with INQuant(predictions_file=predictions_file,
             mzml_file_list=mzml_files,
             proteome_file=proteome_file,
             script_purpose='Mixed proteome quantification with MBR',
             experiment_name='mixed_proteome_MBR',
             ppm_tolerance=new_tolerance,
             confidence_filter=0.0,
             output_file_path=output_folder) as mixed_MBR:
    
    mixed_MBR.run(mbr=True,
            overwrite_all=True)

# Running without MBR
with INQuant(predictions_file=predictions_file,
             mzml_file_list=mzml_files,
             proteome_file=proteome_file,
             script_purpose='Mixed proteome quantification without MBR',
             experiment_name='mixed_proteome_normal',
             ppm_tolerance=new_tolerance,
             confidence_filter=0.0,
             output_file_path=output_folder) as mixed_normal:
    
    mixed_normal.run(mbr=False,
            overwrite_all=True)

