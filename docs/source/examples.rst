=========================
Example Usage
=========================

Here we will provide examples of different ways to use both the INQuant and Calibration classes.
Example data can be genereated with the ``data_generation.py`` script.

The example data generator will generate a predictions file, 3 mzML files named file_1.mzml, file_2.mzml, and file_3.mzml, and a proteome file named example_proteome.fasta.
The data will be generated in the directory, where the script is located, and examples will assume that the inquant_tools module is located in the same directory.

Calibration
==================

At the moment, calibration must be run within the INQuant class, this will be deprecated when the calibration is implemented in the InstaNovo pipeline.
As of now, calibration can be run standalone or combined with the quantification.

Standalone:

.. code-block:: python

    from inquant_tools import INQuant, Calibration  

    with INQuant(predictions_file="predictions.csv",
                     mzml_file_list = ["file_1.mzml", "file_2.mzml", "file_3.mzml"],
                     confidence_filter = 0
                    ) as INQ:          
        # Load predictions using INQuant
        predictions_dict = INQ.load_predictions()

        # Create an instance of Calibration and pass the INQuant predictions
        calibration = Calibration(predictions_dict)

        mzml_files = []
        for mzml_id, mzml_values in INQ.mzml_file_dict.items():
            path = mzml_values['file path']
            name = mzml_values['file name']
            calibration.load_mzml(path)
            calibration.calibration(mzml_id, set_training_seed=42)
            new_predictions, new_mzml = calibration.write_calibrated_file(INQ.experiment_path, mzml_name=mzml_id, predictions_name=INQ.experiment_name)
            mzml_files.append(new_mzml)
        
        new_tolerance = calibration.get_ppm_tolerance()

        print("New tolerance: ", new_tolerance)

The output is the maximum difference of error between the theoritical and calibrated m/z values in ppm. This can be used as the ``ppm_tolerance``, in INQuant. 
The script will also generate new mzML files with the calibrated m/z values, until this is deprecated, when the calibration is implemented in the InstaNovo pipeline.

INQuant
==================

Using the ``run`` function:
--------------------------------
No user arguments:

.. code-block:: python

    from inquant_tools import INQuant  

    with INQuant(predictions_file="predictions.csv",
                 mzml_file_list = ["file_1.mzml", "file_2.mzml", "file_3.mzml"],
                 proteome_file="example_proteome.fasta"
                ) as INQ:          
        INQ.run()

INQuant will promt the user for the column names of the predictions file, if these do not correspond to the default names.
The ``run`` function combines MBR, quantification, protein alignment and inference, and output of files. All of these can be configured with the user arguments or run seperately.


Combined calibration and quantification:
=========================================

.. literalinclude:: examples/example_script.py
   :language: python


