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

proteome_file = path_to_input_files + 'HumanProteome.fasta'

output_folder = 'Results/GluC/'

# Running with MBR
with INQuant(predictions_file=predictions_file,
             mzml_file_list=mzml_files,
             proteome_file=proteome_file,
             script_purpose='GluC quantification with MBR',
             experiment_name='GluC_MBR',
             ppm_tolerance=49.25442393109819,
             confidence_filter=0.0,
             output_file_path=output_folder) as GluC_MBR:
    
    GluC_MBR.run(mbr=True,
            overwrite_all=True)

# Running without MBR
with INQuant(predictions_file=predictions_file,
             mzml_file_list=mzml_files,
             proteome_file=proteome_file,
             script_purpose='GluC quantification without MBR',
             experiment_name='GluC_normal',
             ppm_tolerance=49.25442393109819,
             confidence_filter=0.0,
             output_file_path=output_folder) as GluC_normal:
    
    GluC_normal.run(mbr=False,
            overwrite_all=True)

