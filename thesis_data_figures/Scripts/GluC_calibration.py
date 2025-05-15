from inquant_tools import *

path_to_input_files = 'Datafiles/GluC/'

predictions_file = path_to_input_files + 'gluc_labelled_kpreds.csv'
mzml_files = [
    path_to_input_files + 'mzml/20230802_EV_KK_Evo_Whisper100_20SPD_FAIMS_2CV_1475_Ctrl1.mzML',
    path_to_input_files + 'mzml/20230802_EV_KK_Evo_Whisper100_20SPD_FAIMS_2CV_1475_Ctrl2.mzML',
    path_to_input_files + 'mzml/20230802_EV_KK_Evo_Whisper100_20SPD_FAIMS_2CV_1475_Ctrl3.mzML',
    path_to_input_files + 'mzml/20230802_EV_KK_Evo_Whisper100_20SPD_FAIMS_2CV_1475_GluC1.mzML',
    path_to_input_files + 'mzml/20230802_EV_KK_Evo_Whisper100_20SPD_FAIMS_2CV_1475_GluC2.mzML',
    path_to_input_files + 'mzml/20230802_EV_KK_Evo_Whisper100_20SPD_FAIMS_2CV_1475_GluC3.mzML'
]

output_folder = 'Results/GluC/'

# Making the calibrated mzML files
with INQuant(predictions_file=predictions_file,
             mzml_file_list=mzml_files,
             script_purpose='GluC calibration with labelled data',
             output_file_path=output_folder) as INQ:
    
    # Load predictions using INQuant
    predictions_dict = INQ.load_predictions()

    calibration = Calibration(predictions_dict)

    for mzml_id, mzml_values in INQ.mzml_file_dict.items():
        path = mzml_values['file path']
        name = mzml_values['file name']
        calibration.load_mzml(path)
        calibration.calibration(mzml_id, set_training_seed=42)
        new_mzml = calibration.write_calibrated_mzml_file(output_folder, mzml_name=mzml_id)
        INQ.status_file_creation.append(f"- Wrote calibrated mzML file: {new_mzml}\n")
        mzml_max_tolerance = calibration.ppm_tolerances[mzml_id]
        INQ.variable_status.append(f"- Maximum ppm tolerance for {mzml_id}: {mzml_max_tolerance}\n")

    new_tolerance = calibration.get_ppm_tolerance()
    INQ.variable_status.append(f"- Maximum ppm tolerance across all experiments: {new_tolerance}\n")
   

