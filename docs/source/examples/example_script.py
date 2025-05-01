from inquant_tools import INQuant, Calibration

mzml_file_list=[f"file_{num}.mzML" for num in [1,2,3]]
predictions = "predictions.csv"

with INQuant(
    predictions_file=predictions,
    mzml_file_list=mzml_file_list,
    script_purpose='Calibration example and runthrough of INQuant with auto generated example data',
    confidence_filter=0.0
    ) as INQ:
    # Load predictions using INQuant
    predictions_dict = INQ.load_predictions()

    # Create an instance of Calibration and pass the INQuant instance
    calibration = Calibration(predictions_dict)

    mzml_files = []
    for mzml_id, mzml_values in INQ.mzml_file_dict.items():
        path = mzml_values['file path']
        name = mzml_values['file name']
        calibration.load_mzml(path)
        calibration.calibration(mzml_id, set_training_seed=42)
        new_mzml = calibration.write_calibrated_mzml_file(INQ.experiment_path, mzml_name=mzml_id, predictions_name=INQ.experiment_name)
        mzml_files.append(new_mzml)
    
    new_tolerance = calibration.get_ppm_tolerance()

    INQ.ppm_tolerance = new_tolerance

    INQ.mzml_file_list = mzml_files

    INQ.run(mbr = True, 
            overwrite_all=True)
