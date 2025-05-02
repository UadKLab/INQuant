from pyopenms import *
import pandas as pd
import numpy as np
import time
from Bio import SeqIO
import re
from tqdm.auto import tqdm
# From the module import functions, classes and variables
from . import *
 
class INQuant():

    def __init__(self,
                 predictions_file, 
                 mzml_file_list,
                 proteome_file = None,  
                 script_purpose = None, 
                 ppm_tolerance = 50,
                 empty_values = '',  
                 confidence_filter = 0.95,
                 experiment_name = str(), 
                 output_file_path = str(), 
                 intensity_tolerance = 0.2, 
                 rt_tolerance = 0.01):
        """ Class to quantify and group peptides predicted by InstaNovo.

        It is recommended to always initialize the INQuant class in a 
        with-statement. This will provide a status file upon completion of the 
        script with an overview of time usage and computations completed.

        Parameters
        ----------
        predictions_file : str
            Path to the predictions file.
        mzml_file_list : list
            List of paths to mzML files.
        proteome_file : str, optional
            Path to the proteome file. Optional but needed for alignment and any of the protein 
            grouping functions.
        script_purpose : str, optional
            Purpose of the script.
        ppm_tolerance : float, optional
            Mass and mass charge tolerance for the quantification in ppm. Default is 50.
        empty_values : str, optional
            How empty fields in all output files will be filled. Default is '' (empty), 
            but can be set to '0', 'NA', or any other string.
        confidence_filter : float, optional
            Confidence filter for the predictions. Default is 0.95.
        experiment_name : str, optional
            Name of the experiment. Will be used for file names if specified; 
            otherwise, it uses the predictions filename.
        output_file_path : str, optional
            Path to the output files. Default is None. 
        intensity_tolerance : float, optional
            Relative tolerance for intensity. When intensity falls below the relative tolerance compared to the intensity of the identification scan, it sets the quantification boundary.
            Default is 0.2.
        rt_tolerance : float, optional
            Tolerance for retention time window on either side in percent. Default is 0.01.

        Attributes
        ----------
        exception_called : bool
            Flag indicating if an exception was raised during the script execution.
        break_message : str
            Message associated with the exception raised.
        start : float
            Start time of the script execution.
        predictions_file : str
            Path to the predictions file.
        experiment_name : str
            Name of the experiment.
        experiment_path : str
            Path to the output files.
        status_file : str
            Path to the status file.
        mzml_file_dict : dict
            Dictionary containing mzML file information.
        mzml_file_list : list
            List of paths to mzML files.
        unused_mzml : dict
            Dictionary containing unused mzML files.
        confidence_filter : float
            Confidence filter for the predictions.
        ppm_tolerance : float
            Mass and mass charge tolerance for the quantification in ppm.
        intensity_tolerance : float
            Relative tolerance for intensity.
        rt_tolerance : float
            Tolerance for retention time window on either side in percent.
        psm_dict : dict
            Dictionary containing peptide-spectrum match information.
        psm_dict_index : int
            Index for the peptide-spectrum match dictionary.
        empty_values : str
            How empty fields in all output files will be filled.
        unimod_dict : dict
            Dictionary containing UniMod information to limit multiple lookup.
        status_file_creation : list
            List of status file creation information.
        status_noise : list
            List of noise for each mzML file.
        timer_counts : list
            List of time flags for the script execution.
        files_opened : set
            Set of files opened during the script execution.
        script_purpose : str
            Purpose of the script.
        variable_status : list
            List of variable status information.
        proteome_file : str
            Path to the proteome file.   

        Raises
        ------
        Break
            If no predictions file is provided, the script will raise a Break exception and terminate.

        """

        # Initializing attributes for status file
        self.start = time.perf_counter()
        self.exception_called = False
        self.break_message = None
        self.start = [time.perf_counter(), time.localtime()]
        self.status_file_creation = []
        self.status_noise = []
        self.timer_counts = []
        self.files_opened = set()
        self.script_purpose = script_purpose
        self.variable_status = []
        
        if not predictions_file: # Needed for initializing the rest for clean exit, but will call Break exception later if no predictions file
            predictions_file = ''

        self.predictions_file = predictions_file
        
        # Initializing names and paths for output
        if experiment_name == '':
            self.experiment_name = self.predictions_file.split('/')[-1].split('.')[0]
        else:
            self.experiment_name = experiment_name
        
        self.experiment_path = output_file_path
        
        if confidence_filter != 0.0:
            self.experiment_name += f"_filtered_{confidence_filter}"
        else:
            self.experiment_name += '_unfiltered'

        self.status_file = f"{self.experiment_path}{self.experiment_name}_status.txt"

        # For load_predictions
        self.mzml_file_dict = {}
        self.mzml_file_list = mzml_file_list
        self.unused_mzml = {}
        self.confidence_filter = confidence_filter
        
        # For make_psm_table
        self.ppm_tolerance = ppm_tolerance

        # For quant_calc
        self.intensity_tolerance = intensity_tolerance
        self.rt_tolerance = rt_tolerance
        self.psm_dict = {}
        self.psm_dict_index = 0
        self.empty_values = empty_values
        self.unimod_dict = {}

        # Protein table
        self.proteome_file = proteome_file

        if predictions_file == '':
            self.__exit__(Break,'Error: No predictions file was provided.',None)
        
        if len(mzml_file_list) == 0 or not mzml_file_list:
            self.__exit__(Break,'Error: No mzml files were provided.',None)


    def __enter__(self):
        """ Prepares the object for use in a with-statement. This method ensures 
        compatibility with the `with` statement and is called when the 
        `with` block is entered.

        Returns
        -------
        self : object
            The object itself, allowing it to be used within the `with` block.
        
        """

        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        """
        Finalizes the status file for the experiment when the `with` block is exited. 
        This method writes the experiment details, including the script's runtime, 
        any encountered errors, and other relevant information to the status file.

        Parameters
        ----------
        exc_type : type or None
            The exception type raised, or `None` if no exception was raised.

        exc_value : Exception or None
            The exception instance raised, or `None` if no exception was raised.

        traceback : traceback or None
            The traceback object, or `None` if no exception was raised.

        Attributes
        ----------
        stop : float
            The end time of the script execution.

        Returns
        -------
        None
            This method does not return any value. It writes details to the status file.
        """
        
        # Checks reason for exit
        if exc_type is Break or self.exception_called:
            if exc_value:
                print(exc_value)
                if self.break_message == None: # Only the first error message will be written to the status file
                    self.break_message = exc_value 
            self.exception_called = True
        
        # If no exception was raised, write the status file
        if exc_type is None: 
            self.stop = time.perf_counter()
            time_elapsed = self.stop - self.start[0]
            if time_elapsed > 180:
                time_string = f"{time_elapsed/60:.3} minutes"
            else:
                time_string = f"{time_elapsed:.4} seconds"
        
            status = open(self.status_file,'w')
            
            # Headers
            status.write(f"Status for script run for experiment: {self.experiment_name}\n\n")
            if self.script_purpose != '':
                status.write(f"Script purpose:\n{self.script_purpose}\n\n")
            if self.break_message != None:
                status.write(f"Script broke during execution due to: \n{self.break_message}\n\n")
            
            # Add date and time of script run when script started
            status.write(f"Script started at: {time.strftime('%Y-%m-%d %H:%M:%S',self.start[1])}\n\n")

            status.write(f"Time spent on entire script: {time_string}\n")
            if self.timer_counts != []:
                status.write(f"- Time elapsed for time-flags called underway in the script:\n")
                for flag in self.timer_counts:
                    status.write(flag)

            # Add all variables for calculations
            status.write(f"\nMass charge tolerance: {self.ppm_tolerance} ppm\n")
            status.write(f"Confidence filter: {self.confidence_filter*100}%\n") if self.confidence_filter != 0.0 else status.write(f"Confidence filter: None\n")
            status.write(f"Intensity tolerance: {self.intensity_tolerance*100}%\n")
            status.write(f"Retention time window tolerance: {self.rt_tolerance*100}%\n")
            if len(self.files_opened) != 0:
                status.write(f"Empty values in output files: {self.empty_values}\n")
            if self.variable_status != []:
                for elm in self.variable_status:
                    status.write(f"{elm}")
            
            # All files opened or written in the script
            status.write('\nFiles loaded in this script (excluding mzML):\n')

            status.write(f"-- Predictions file: {self.predictions_file}\n")
            with open(self.predictions_file, "r") as f:
                pred_lines = sum(1 for _ in f)
            status.write(f"---- Line count: {pred_lines}\n")
            
            if len(self.files_opened) != 0:
                for line in self.files_opened:
                    status.write(line)

            status.write('\nmzML files loaded in running of this script:\n')
            if self.mzml_file_dict != {}:
                for key,val in self.mzml_file_dict.items():
                    status.write(f"ID: {key}, filename: {val['file name']}\n")
            else:
                status.write('None were loaded.\n')
            
            # Noise for each mzML file 
            if self.status_noise != []:
                status.write(f"\nNoise calculations in this script based on {self.noise_variables}\n")
                for elm in self.status_noise:
                    status.write(f"{elm}\n")
            
            if self.unused_mzml != {}: 
                status.write('\nmzML files not loaded to this script from same dataset:\n')
                for key in self.unused_mzml.keys():
                    status.write(f"ID: {key}, name in sample column: {self.unused_mzml[key]}\n")

            if self.status_file_creation != []:
                status.write(f"\nIn this script, the following data files were created:\n")
                for elm in self.status_file_creation:
                    status.write(elm)

            status.close()
            print(f"The experiment has run, details can be found in {self.status_file}. Remember to rename the status file if you would like to keep it before running next script in this experiment.")

    def INQuant_close(self):
        """
        Exit function writes the status file for the experiment, should be run if the script is run without using the with-statement. 
        """
        self.__exit__(None,None,None)

    def time_flag(self, name):
        """ Records the time taken since the start of the script and appends the 
        result to the timer counts. If a `time_flag` is set anywhere within 
        the with-statement, the status file will show how long the script took 
        to reach the `time_flag`.

        Parameters
        ----------
        name : str
            The name associated with the time flag. This will be used in the status 
            file to indicate the specific point in the script.

        Returns
        -------
        None
            This method does not return any value. It appends the time information 
            to the `timer_counts` list attribute.
        
        """

        time_stop = time.perf_counter()
        time_elapsed = time_stop - self.start[0]
        if time_elapsed > 180:
            time_string = f"{time_elapsed/60:.3} minute(s)\n"
        else:
            time_string = f"{time_elapsed:.4} second(s)\n"

        self.timer_counts.append(f"-- {name}: {time_string}")

    def run(self,  
            noise_boundary = 0.05,
            noise_ms_level = 2,
            mbr = False,
            mbr_tolerance = 0.9,
            cleavage_length = 4,
            quant_method = 'mean', 
            top_n_peptides = 5, 
            normalize_abundance = 'median',
            write_psm_table = True, 
            psm_file_name = str(), 
            write_peptide_table = True, 
            peptide_file_name = str(),
            write_protein_table = True, 
            protein_file_name = str(),
            write_ungrouped_protein_table = False, 
            ungrouped_protein_file_name = str(),
            write_individual_quantification_files = False,
            individual_quant_file_name = str(),
            overwrite_all = False):
        
        """ Perform the full quantification process by combining multiple functionalities 
        of the class into a single command. This should be executed when the class is 
        initialized with a `with` statement to ensure proper management of resources.

        Parameters
        ----------
        noise_boundary : float, optional
            The boundary for noise calculation. Default is 0.05.
        noise_ms_level : int, optional
            The MS level for noise calculation. Default is 2.
        mbr : bool, optional
            Whether to run the MBR (Match Between Runs) algorithm. Default is False.
        mbr_tolerance : float
            The minimum confidence for the MBR spectras. Only predictions with a higher confidence
            will have mbr matches. Default is 0.9.
        cleavage_length : int, default=4
            The length of the cleavage site used for alignment. This parameter is used to output the protein position, 
            with cleavage length being the number of amino acids before and after the peptide sequence in the protein.
        quant_method : str, optional
            The quantification method to use. Default is 'mean'.
        top_n_peptides : int, optional
            The number of top peptides to consider for quantification values for each protein. Default is 5.
        normalize_abundance : str, optional
            The method for normalizing abundance values. Options include 'median', 
            'mean', 'tic', or 'false'. Default is 'median'.
        write_protein_table : bool, optional
            Whether to write the protein table to a file. Default is True.
        protein_file_name : str, optional
            The name of the protein table file. Default is '[experiment_name]_protein_table.csv'.
        write_peptide_table : bool, optional
            Whether to write the peptide table to a file. Default is True.
        peptide_file_name : str, optional
            The name of the peptide table file. Default is '[experiment_name]_peptide_table.csv'.
        write_psm_table : bool, optional
            Whether to write the PSM (Peptide-Spectrum Match) table to a file. Default is True.
        psm_file_name : str, optional
            The name of the PSM table file. Default is '[experiment_name]_psm_table.csv'.
        write_ungrouped_protein_table : bool, optional
            Whether to write the ungrouped protein table to a file. Default is False.
        ungrouped_protein_file_name : str, optional
            The name of the ungrouped protein table file. Default is '[experiment_name]_ungrouped_protein_table.csv'.
        write_individual_quantification_files : bool, optional
            Whether to write individual quantification files for each mzML file. Default is False.
        individual_quant_file_name : str, optional
            The name of the individual quantification file. Default is '[experiment_name]_quantification_[mzml_file_id].csv'.
        overwrite_all : bool, optional
            Whether to overwrite all existing files without prompting. Default is False, 
            which will ask for confirmation before overwriting any file.

        Raises
        ------
        ValueError
            If the `top_n_peptides` parameter is less than 1, a `ValueError` will be raised.
        
        """
        
        # If protein quantification is to be performed, break if no peptides can be used for quantification
        if top_n_peptides < 1 and self.proteome_file != None:
            self.__exit__(Break,"ValueError: Top N peptides must be greater than 0. No files were written.",None)
        
        # Always check if exception was called, so class will gracefull exit if error was raised
        if self.exception_called:
            return
        
        # Load predictions with MBR specifications
        predictions_dict = self.load_predictions(mbr=mbr, mbr_tolerance=mbr_tolerance)
        
        # Make the PSM table 
        psm_dict = self.make_psm_table(predictions_dict=predictions_dict, noise_boundary=noise_boundary, noise_ms_level=noise_ms_level)
        
        # Write the individual quantification files if specified
        self.write_files(write_protein_table=False, 
                        write_peptide_table=False,
                        write_psm_table=False,
                        write_individual_quantification_files = write_individual_quantification_files,                            
                        individual_quant_file_name = individual_quant_file_name,
                        overwrite_all = overwrite_all)
        
        # Make the peptide table from the PSM table 
        peptide_dict = self.make_peptide_table(psm_dict)

        # If there are multiple mzml files and normalization is specified, normalize the abundance values in the peptide table 
        if len(self.mzml_file_list) > 1 and normalize_abundance.lower() != 'false':
            self.normalizer(peptide_dict, type = normalize_abundance)

        # Compute protein alignment, quantification and inference if a proteome file is provided
        if self.proteome_file != None:
            self.load_proteome(self.proteome_file)
            self.compute_alignments(cleavage_length=cleavage_length)
            self.make_protein_table(top_n_peptides=top_n_peptides, quant_method=quant_method)
            self.write_files(write_protein_table = write_protein_table,
                            protein_file_name = protein_file_name,
                            write_peptide_table = write_peptide_table,
                            peptide_file_name = peptide_file_name,
                            write_psm_table = write_psm_table,
                            psm_file_name = psm_file_name,
                            write_ungrouped_protein_table = write_ungrouped_protein_table,
                            ungrouped_protein_file_name = ungrouped_protein_file_name,
                            overwrite_all = overwrite_all)
            
        else:
            # If no proteome file is provided, write the peptide table and PSM table only
            self.write_files(write_protein_table = False,
                            write_peptide_table = write_peptide_table,
                            peptide_file_name = peptide_file_name,
                            write_psm_table = write_psm_table,
                            psm_file_name = psm_file_name,
                            overwrite_all = overwrite_all)


    def load_predictions(self,
                         mbr = False, 
                         mbr_tolerance = 0.9) -> dict:
        """ Loads predictions into the class and assigns file IDs for all mzml files in the experiment. 
        If only one mzml file is provided, the file ID will always be 'mzml'. This function handles 
        multiple mzml files and attempts to match their identifiers to those in the predictions file. 
        Additionally, if mbr (Match Between Runs) is enabled, it processes missing data by 
        imputing values from previous replicates.

        Parameters
        ----------
        mbr : bool, optional, default is False
            Flag indicating whether to run mbr (Match Between Runs). If True, missing values 
            for certain replicates will be imputed based on available data from other replicates.

        mbr_tolerance : float, default is 0.9
            The minimum confidence for the mbr specs. Only predictions with a higher confidence
            will have mbr matches. 

        Returns
        -------
        sorted_predictions_dict : dict 
            The function modifies internal class attributes such as `mzml_file_dict`, `sorted_predictions_df`,
            and `sorted_predictions_dict`. If `mbr` is True, the method also populates the mbr specifications 
            into `sorted_predictions_dict`.
            
        Raises
        ------
        Break : 
            If there is an error matching mzml files with prediction IDs or other issues arise, the process 
            will stop and the Break exception will be raised.
        
        """
        # Always check if exception was called, so class will gracefull exit if error was raised
        if self.exception_called:
            return

        # Ensure multiple mzML files before allowing MBR
        if len(self.mzml_file_list) <= 1:
            mbr = False
            print('There was only one mzML file provided. MBR will not be run.')
        

        def find_mbr(df = pd.DataFrame, mbr_tolerance = 0.9):
            
            # Load identifiers for the different mzML files
            mzml_dict = self.mzml_file_dict

            # Initialize columns for spectra data
            df['mz'] = np.nan
            df['rt'] = np.nan
            df['charge'] = np.nan
            # Initialize dict for storing experiment lengths for relative RT adjustments
            exp_len = {}

            # Populate rt and mz from mzml files
            for _, mzml_id in tqdm(enumerate(mzml_dict.keys()), 
                                   total=len(mzml_dict.keys()), 
                                   desc='Loading mzML files for MBR'):
                # For each mzML file, load the spectra data
                mask = df['id'] == mzml_id
                exp = MSExperiment()
                MzMLFile().load(mzml_dict[mzml_id]['file path'], exp)
                spectra_dict = {
                    i: {
                        'mz': exp.getSpectrum(i).getPrecursors()[0].getMZ() if exp.getSpectrum(i).getPrecursors() else np.nan,
                        'rt': exp.getSpectrum(i).getRT(),
                        'charge': exp.getSpectrum(i).getPrecursors()[0].getCharge() if exp.getSpectrum(i).getPrecursors() else ''
                    }
                    for i in range(exp.getNrSpectra())
                }
                # Store the highest RT (experiment length)
                exp_len[mzml_id] = spectra_dict[exp.getNrSpectra() - 1]['rt']

                # Add the spectra data to the predictions file
                for i in df[mask].index:
                    scan_number = df.at[i, 'scan_number']
                    df.at[i, 'rt'] = spectra_dict[scan_number]['rt']
                    df.at[i, 'mz'] = spectra_dict[scan_number]['mz']
                    df.at[i, 'charge'] = spectra_dict[scan_number]['charge']
                    df.at[i, 'mass'] = df.at[i, 'mz'] * df.at[i, 'charge'] - (proton_da * df.at[i, 'charge'])

            # Compare across all files (not by replicate)
            complete_df_rows = []
            fully_identified_sequences = []
            all_ids = list(mzml_dict.keys())
            df['group'] = 0
            
            # Group by 'preds' and mass within +- ppm_tolerance, assigning group numbers for easier iteration
            sequence_groups = df.groupby(by = ['preds'])
            group_number = 1
            for sequence, group in sequence_groups:
                group.sort_values(by=['mass'], ascending=True, inplace=True)
                prev_mass = 0
                for i, row in group.iterrows():
                    if prev_mass == 0:
                        df.at[i, 'group'] = group_number
                    else:
                        if abs(row['mass'] - prev_mass) <= (self.ppm_tolerance/10**6)*prev_mass:
                            df.at[i, 'group'] = group_number
                        else:
                            group_number += 1
                            df.at[i, 'group'] = group_number

                    prev_mass = row['mass']
                
                group_number += 1 

            sequence_groups = df.groupby(by = ['group'])
            
            # Iterate through groups and assign MBR spectras where missing values
            for _, (group_no, group) in tqdm(enumerate(sequence_groups), 
                                             total=len(sequence_groups), 
                                             desc='Retrieving MBR spectras'):
                # Select best row per file (highest log_prob)
                complete_df_rows.append(group.copy())
                # Sort based on log_probs to use the best row as the reference row
                group.sort_values('log_probs', ascending=False, inplace = True)
                # Drop duplicates if there are multiple identifications of the same peptide within tolerance in the group
                # Best confidence will be kept
                group.drop_duplicates(subset='id', inplace = True)
                found_ids = group['id'].tolist()
                # Find the mzML IDs with missing identifications
                missing_ids = [f_id for f_id in all_ids if f_id not in found_ids]

                if missing_ids:
                    # There must be one in the group with a score above the tolerance
                    if not any(group['log_probs'] > np.log10(mbr_tolerance)):
                        continue
                    
                    # Find previous sequence with complete data to use as reference 
                    try:
                        prior_group = fully_identified_sequences[-1]
                    except IndexError:
                        continue
                    
                    # Calculate rt deltas
                    deltas = []
                    for f_id in found_ids:
                        try:
                            curr_rt = group[group['id'] == f_id].iloc[0]['rt']
                            prev_rt = prior_group[prior_group['id'] == f_id].iloc[0]['rt']

                            deltas.append((curr_rt - prev_rt)/exp_len[prior_group['id'].iloc[0]])
                        except IndexError:
                            continue
                    avg_delta = np.mean(deltas) if deltas else 0.0

                    # Base mass charge and charge on the row with best confidence in the group
                    best_row = group[group['log_probs'] == max(group['log_probs'])].iloc[0]
                    
                    for m_id in missing_ids:
                        try:
                            prev_row = prior_group[prior_group['id'] == m_id].iloc[0].copy()
                            new_row = prev_row.copy()
                            new_row['preds'] = group['preds'].iloc[0]
                            new_row['mz'] = best_row['mz']
                            new_row['scan_number'] = None
                            new_row['charge'] = best_row['charge']
                            new_row['rt'] = prev_row['rt'] + avg_delta*exp_len[m_id]
                            new_row['log_probs'] = None
                            new_row['group'] = group['group'].iloc[0]
                            complete_df_rows.append(pd.DataFrame([new_row]))
                        except IndexError:
                            continue
                else: # All mzML files have a identification for this peptide
                    fully_identified_sequences.append(group)

            # Combine and return all rows
            final_df = pd.concat(complete_df_rows, ignore_index=True)
            final_df = final_df.sort_values(by=['group']).reset_index(drop=True)

            # Convert to dictionary with dataframes for each mzML file to be used in quantification
            final_dict = {
                mzml_id: final_df[final_df['id'] == mzml_id].drop(columns=['id', 'group']).reset_index(drop=True)
                for mzml_id in mzml_dict.keys()
            }

            self.time_flag(f"Finished running MBR")
            
            return final_dict
        
        # Finding the shortest possible identifiers for the mzML files and assign these as IDs
        def find_ids(list_of_strings):
            names = [str(string).split('/')[-1].split('.')[0] for string in list_of_strings]
            min_name = min([len(name) for name in names])

            unique_list_slice_start = 0
            unique_list_slice_stop = 0

            for i in range(min_name):
                unique_char = set([string[i] for string in names])
                if len(unique_char) != 1:
                    unique_list_slice_start = i
                    break

            for i in range(1, min_name+1):
                i *= -1
                unique_char = set([string[i] for string in names])
                if len(unique_char) != 1:
                    unique_list_slice_stop = i+1
                    break
            
            # Make a final list of the unique IDs
            unique_id = [string[unique_list_slice_start:len(string)+unique_list_slice_stop] for string in names]
            
            return unique_id
        
        sorted_predictions_df = pd.read_csv(self.predictions_file, sep=',')
        sorted_predictions_dict = {}
        
        # Check the column names in the predictions file
        expected_columns = ['scan_number', 'preds', 'log_probs', 'file']
        while True:
            missing_columns = [col for col in expected_columns if col not in sorted_predictions_df.columns]
            if missing_columns:
                print(f"Error: Missing columns in predictions file: {', '.join(missing_columns)}")
                answ = input("Do you want to provide the column names in the predictions file? (y/n): ").strip().lower()
                if answ == 'y':
                    for col in missing_columns:
                        col_name = input(f"Please provide the column name in the predictions file, which corresponds to '{col}': ")

                        sorted_predictions_df.rename(columns={col_name: col}, inplace=True)
                else:
                    print("Exiting the program.")
                    self.__exit__(Break,'Error: Missing columns in predictions file.',None)
                    return
            else:
                break


        # Filtering out predictions with low confidence
        sorted_predictions_df = sorted_predictions_df[np.power(10,sorted_predictions_df['log_probs'].astype(float)) >= self.confidence_filter]
        sorted_predictions_df.dropna(subset=['preds'], inplace=True)
        
        sorted_predictions_df['scan_number'] -= 1  # To match with 0-indexing from mzML files

        try: 
            samples = sorted_predictions_df['file'].drop_duplicates().to_list() # Finding unique sample names/experiments in predictions file
            unique_predictions_id = find_ids(samples)
            if len(self.mzml_file_list) == 1: 
                mzml_file_name = self.mzml_file_list[0].split('/')[-1].split('.')[0]
                possible_ids = [mzml_id if mzml_id in mzml_file_name else '' for mzml_id in unique_predictions_id]

                if(sum(bool(mzml_id) for mzml_id in possible_ids) == 0):
                    print(f"IDs found in predictions file do not match mzml-filenames. IDs found in predictions file:\n{unique_predictions_id}\nmzML-filename:\n{mzml_file_name}")
                    self.__exit__(Break,'Error: IDs found in predictions file do not match mzml-filenames. ',None)
                    return
                
                else:
                    max_length = max(len(mzml_id) for mzml_id in possible_ids)
                    longest_ids = [mzml_id for mzml_id in possible_ids if len(mzml_id) == max_length]

                    # Check if there are multiple ids with the same maximum length
                    if len(longest_ids) > 1:
                        print(f"IDs found in predictions file match multiple mzml-filenames.\n{unique_predictions_id}\nIdetified following possible ID matches:\n{longest_ids}\nmzML file name:\n{mzml_file_name}")
                        self.__exit__(Break,'Error: IDs found in predictions file do not match mzml-filenames. ',None)
                        return

                    pred_index = unique_predictions_id.index(longest_ids[0])
                    mzml = longest_ids[0]
                    
                    mask = sorted_predictions_df['file'] == samples[pred_index]
                    sorted_predictions_df.loc[mask, 'id'] = mzml # Adding ID to the predictions file

                    self.mzml_file_dict[mzml] = {} # Making dict for mzml file name types
                    self.mzml_file_dict[mzml]['file path'] = self.mzml_file_list[0] # Adding file path to mzml file
                    self.mzml_file_dict[mzml]['file name'] = mzml_file_name # Adding file name to mzml file
                    
                    # Remove the used mzml file and prediction ID from the lists
                    unique_predictions_id.pop(pred_index)
                    samples.pop(pred_index)

            else: 
                unique_mzml_id = find_ids(self.mzml_file_list)
                for (mzml_index, mzml) in enumerate(unique_mzml_id):
                    try: 
                        pred_index = unique_predictions_id.index(mzml)
                    except ValueError:
                        print(f"IDs found in predictions file do not match mzml-filenames. IDs found in predictions file:\n{unique_predictions_id}\nID found in mzml-filename:\n{unique_mzml_id}")
                        self.__exit__(Break,'Error: IDs found in predictions file do not match mzml-filenames. ',None)
                        return
                    
                    mask = sorted_predictions_df['file'] == samples[pred_index]
                    sorted_predictions_df.loc[mask, 'id'] = mzml # Adding ID to the predictions file
                                                                
                    self.mzml_file_dict[mzml] = {} # Making dict for mzml file name types
                    self.mzml_file_dict[mzml]['file path'] = self.mzml_file_list[mzml_index] # Adding file path to mzml file
                    self.mzml_file_dict[mzml]['file name'] = self.mzml_file_list[mzml_index].split('/')[-1].split('.')[0] # Adding file name to mzml file
                    
                    # Remove the used mzml file and prediction ID from the lists
                    unique_predictions_id.pop(pred_index)
                    samples.pop(pred_index)

            if(len(unique_predictions_id) != 0):
                for mzml_id in unique_predictions_id:
                    self.unused_mzml[mzml_id] = samples[unique_predictions_id.index(mzml_id)]
        
        except (KeyError, IndexError):
            filtered_sorted_predictions_df = sorted_predictions_df.copy()
            sorted_predictions_dict['mzml'] = filtered_sorted_predictions_df
            
            self.mzml_file_dict['mzml'] = {} # Making dict for mzml file name types
            self.mzml_file_dict['mzml']['file path'] = self.mzml_file_list[0] # Adding file path to mzml file
            self.mzml_file_dict['mzml']['file name'] = self.mzml_file_list[0].split('/')[-1].split('.')[0] # Adding file name to mzml file

        # Load the predictions into the dictionary, with or without MBR
        self.variable_status.append(f"Ran MBR: {mbr}\n")
        
        if mbr == True:
            self.variable_status.append(f"MBR tolerance: {mbr_tolerance*100}%\n")
            sorted_predictions_dict = find_mbr(sorted_predictions_df, mbr_tolerance=mbr_tolerance)
        
        else:
            for mzml_id in self.mzml_file_dict.keys():
                mask = sorted_predictions_df['id'] == mzml_id
                sorted_predictions_dict[mzml_id] = sorted_predictions_df[mask].drop(columns=['id']).reset_index(drop=True, inplace = False)

        return sorted_predictions_dict
            
    def make_psm_table(self, 
                       predictions_dict = dict(), 
                       noise_boundary = 0.05, 
                       noise_ms_level = 2) -> dict:
        """ Creates a Peptide-Spectrum Match (PSM) table from the predictions data.
        This function is dependent on the `load_predictions` method and requires the
        predictions data to be loaded first, because of the mzML IDs which are generated
        in the `load_predictions` method.
        Quantifications use user specifications from the initialization of the class.

        Parameters
        ----------
        predictions_dict : dict
            Dictionary containing the predictions data. The keys are mzml file IDs and the 
            values are DataFrames containing the predictions for each file.
        
        noise_boundary : float, optional, default is 0.05
            Boundary value for noise calculation. This value will be used to define the 
            threshold below which signals are considered as noise.

        noise_ms_level : int, optional, default is 2
            MS level to be used for noise calculation. Determines which level of the 
            mass spectrometry data (MS1 or MS2) will be considered for noise analysis.

        Attributes
        ----------
        noise_variables : str
            String containing the noise calculation parameters used in the experiment.
        
        Returns
        -------
        psm_dict : dict
            Dictionary containing the PSM data with the calculated quantification values.

        Notes
        -----
            - The noise calculation is based on the provided boundary and MS level, which helps in filtering out irrelevant signals.
            - Noise for each experiment can be found in the status file, if such a file is generated by using the with-statement for the class.

        """
        # Always check if exception was called, so class will gracefull exit if error was raised
        if self.exception_called:
            return
        
        # Function for noise calculation for each experiment
        def get_noise_for_whole_exp(specs = dict, noise_boundary = 0.05, noise_ms_level = 2):

            if noise_boundary == 0 or noise_boundary == None:
                return 0
            
            spec_data = [peak for elm in specs.values() if elm['ms_level'] == noise_ms_level for peak in elm['peaks'][1]]

            # Find lowest 5 percent:
            total_length = len(spec_data)

            amount_for_noise = int(total_length*noise_boundary)
            spec_data.sort()

            spectra_minimums = spec_data[0:amount_for_noise]
            noise = np.mean(spectra_minimums)

            return noise
      
        # For status file
        self.noise_variables = f"MS{noise_ms_level} lowest {noise_boundary*100}% intensities:"
        
        # Function to calculate the quantification values for each mzML file
        def get_quants(spectra_dict = dict, mzml_id = str, noise = float):

            def area_calc(x_axis, y_axis):
                # Calculate the area under the curve
                dx = np.diff(x_axis)
                y_min = np.minimum(y_axis[:-1], y_axis[1:])
                y_diff = np.abs(np.diff(y_axis))
                square_area = dx * y_min
                triangle_area = 0.5 * dx * y_diff

                return np.sum(square_area + triangle_area)
            
            def boundaries(intensities, times, rt_detection_scan):
                # Assign boundaries for the quantification window using the initially specified intensity tolerance 
                
                if len(times) <= 1:
                    return

                if not times[0] < rt_detection_scan < times[-1]:
                    return
                
                high_detect = None

                # Assign preliminary boundaries based on the identification scan
                for i in range(len(times)):
                    if times[i] > rt_detection_scan:
                        low_detect = i-1
                        high_detect = i
                        break
                
                # Ensure large enough window for quantification
                if high_detect is None:
                    return
                
                # Boundaries are set when the intensity falls below the tolerance for two consecutive points
                # The iteration works from the identification scan and outwards in both directions
                id_intensities = intensities[low_detect]
                
                low_boundary = 0
                one = False
                for i in range(low_detect, -1, -1):
                    if intensities[i] < id_intensities * self.intensity_tolerance:
                        if one: 
                            low_boundary = i
                            break
                        else:
                            one = True
                    else: 
                        one = False
    
                high_boundary = len(intensities) - 1
                one = False
                for i in range(high_detect, high_boundary + 1):
                    if intensities[i] < id_intensities * self.intensity_tolerance:
                        if one: 
                            high_boundary = i
                            break
                        else:
                            one = True
                    else:
                        one = False
                
                if not 0 <= low_boundary < high_boundary < len(intensities):
                    return
                
                return low_boundary, high_boundary

            def iterate_spectras(sorted_predictions_dict = dict, noise = float):
                # Main quantification calculation is done here, iteration through each prediction

                filtered_specs = {key: spec for key, spec in spectra_dict.items() if spec['ms_level'] == 1}
                
                for _,  prediction in tqdm(enumerate(sorted_predictions_dict.values()),total=len(sorted_predictions_dict),desc=f"Quantifying peptides for mzml file: {mzml_id}"):
                    
                    # Assign the target mz and RT for the identification
                    target_mz = prediction['meas_mz']
                    rt_first_scan = prediction['rt']
                    
                    # Check if identification is MBR or InstaNovo and assign initial tolerance
                    if prediction['scan_number'] == None:
                        # MBR RT identifications are re-assigned based on the best peak in a 2pct window
                        tolerance_rt = ((1+self.rt_tolerance)*1.02)-1 * rt_first_scan  # Tolerance is large enough to find the best peak and 1 pct. around the best peak in a 2 pct. window
                    else:
                        tolerance_rt = self.rt_tolerance * rt_first_scan
                    
                    # Calculate the specific ppm tolerance
                    tolerance = (target_mz*self.ppm_tolerance)/10**6
                    
                    # Filter the spectra data to only include points within respective tolerances
                    target_time, target_intensities = [], []
                    for index, spec in filtered_specs.items():
                        rt = spec['rt']
                        if abs(rt - rt_first_scan) > tolerance_rt:
                            continue  

                        mz_values = np.array(spec['peaks'][0])
                        Intensity_values = np.array(spec['peaks'][1])
                        
                        mz_mask = np.abs(mz_values - target_mz) < tolerance
                        if len(Intensity_values[mz_mask]) != 0:
                            intensity_sum = np.sum(Intensity_values[mz_mask])
                        else:
                            intensity_sum = None

                        if intensity_sum != None:
                            target_intensities.append(intensity_sum)
                            target_time.append(rt)

                    # Subtract the noise from the intensities
                    target_intensities = [elem - noise for elem in target_intensities]
                    target_intensities = [max(0, elem) for elem in target_intensities]

                    # Convert to numpy arrays for faster calculations
                    target_time = np.array(target_time)
                    target_intensities = np.array(target_intensities)
                    
                    # Only do further calculations if there are relevant points 
                    if target_time.size != 0:
                        # Find the best window for the MBR specs (1pct RT around the max intensity)
                        if prediction['scan_number'] == None:
                            # Finding the max peak in the 2pct window
                            mask = np.abs(target_time - rt_first_scan) <= 2/100 * rt_first_scan
                            target_intensities_2 = target_intensities[mask]
                            if len(target_intensities_2) != 0: # Check if there are any points in the 2pct window
                                target_peak = max(target_intensities_2)
                                target_peak_index = np.where(target_intensities == target_peak)[0][0]
                                # New RT detection scan, find 1pct window around it
                                rt_identify_peak = target_time[target_peak_index]
                                final_mask = np.abs(target_time - rt_identify_peak) <= self.rt_tolerance * rt_identify_peak
                                # Modify the target time and intensities to only include the 1pct window
                                target_time = target_time[final_mask]
                                target_intensities = target_intensities[final_mask]
                                rt_first_scan = rt_identify_peak
                            else: # The RT window is empty, so spec_quant will be NaN when trying to calculate it
                                target_time = target_time[mask]
                                target_intensities = target_intensities_2
                        
                        # Calculate the boundaries and quantification value
                        try:
                            low, high = boundaries(target_intensities, target_time, rt_first_scan)
                            target_time = target_time[low:high + 1]
                            target_ints = target_intensities[low:high + 1]
                            spec_quant = area_calc(target_time, target_ints)
                        
                        # If no boundaries are found or there are no points within the boundaries, error is raised
                        except (TypeError, ValueError):
                            # If the identification was done by MBR, the row will not be inclued
                            if prediction['scan_number'] == None:
                                continue
                            # Else, the quantification will be empty
                            spec_quant = np.nan
                    
                    else:
                        # If no points are found, the quantification will be empty
                        spec_quant = np.nan
                    
                    # Add the prediction to the PSM table
                    self.psm_dict[self.psm_dict_index] = {'sequence': f"{re.sub(r'[^A-Z]', '', prediction['preds'])}".replace('I','L'),
                                                          'seq_modifications': prediction['preds'],
                                                          'no_psms': prediction['no_psms'],
                                                          'mz': prediction['meas_mz'],
                                                          'charge': prediction['charge'],
                                                          'meas_mass': prediction['meas_mass'],
                                                          'rt': rt_first_scan,
                                                          'id_scan': prediction['scan_number'],
                                                          'file_id': mzml_id,
                                                          'exp_file': self.mzml_file_dict[mzml_id]['file name'],
                                                          'conf': prediction['log_probs'], 
                                                          'abundance': spec_quant, 
                                                          'calc_mass': aa_mass['H2O'] + proton_da # H2O, because the AA weights are for bound AA, not free
                                                          }
                    
                    # Add the mass of the sequence to the PSM table
                    for aa in self.psm_dict[self.psm_dict_index]['sequence']:
                        self.psm_dict[self.psm_dict_index]['calc_mass'] += aa_mass[aa]
                    
                    # Add modifications to the PSM table
                    if '(' in self.psm_dict[self.psm_dict_index]['seq_modifications']:
                        # Find all matches along with their starting index
                        mod = [f"{m.start()}({m.group(1)})" for m in re.finditer(r"\((.*?)\)", self.psm_dict[self.psm_dict_index]['seq_modifications'])]
                        
                        for m in mod: 
                            m = re.search(r"\((.*?)\)", m).group(1)
                            try:
                                self.psm_dict[self.psm_dict_index]['calc_mass'] += float(m)
                            except ValueError:
                                if m == 'ox':  
                                    self.psm_dict[self.psm_dict_index]['calc_mass'] += 15.99491463
                                else:
                                    print(f"Error: Could not convert modification mass to float: {m}")

                        self.psm_dict[self.psm_dict_index]['modifications'] = ';'.join(mod)

                    # Implementation for UNIMOD modifications
                    if '[' in self.psm_dict[self.psm_dict_index]['seq_modifications']:
                        # Find all matches along with their starting index
                        mod = [f"{m.start()}[{m.group(1)}]" for m in re.finditer(r"\[(.*?)\]", self.psm_dict[self.psm_dict_index]['seq_modifications'])]
                        
                        for m in mod: 
                            m = re.search(r"\[(.*?)\]", m).group(1)
                            try:
                                self.psm_dict[self.psm_dict_index]['calc_mass'] += self.unimod_dict[m]
                            except KeyError:
                                self.unimod_dict[m] = get_monoisotopic_mass(m)
                                if self.unimod_dict[m] != None:
                                    self.psm_dict[self.psm_dict_index]['calc_mass'] += self.unimod_dict[m]
                                else:
                                    raise TypeError
                            except ValueError or TypeError:
                                print(f"Error: Could not convert modification mass to float: {m}")

                        self.psm_dict[self.psm_dict_index]['modifications'] = ';'.join(mod)
                    
                    # Add a theoretical mz to the PSM table and calculate the error
                    self.psm_dict[self.psm_dict_index]['calc_mz'] = (self.psm_dict[self.psm_dict_index]['calc_mass']+proton_da*(self.psm_dict[self.psm_dict_index]['charge']-1))/self.psm_dict[self.psm_dict_index]['charge']
                    self.psm_dict[self.psm_dict_index]['mz_error'] = self.psm_dict[self.psm_dict_index]['mz'] - self.psm_dict[self.psm_dict_index]['calc_mz']
                    
                    # Add index, so that the PSM table is updated throughout multiple mzML files
                    self.psm_dict_index +=1

            # Filter the data based on the current mzML file ID
            filtered_df = predictions_dict[mzml_id]
            # Add necessary columns to dataframe for quantification
            if 'charge' not in filtered_df.columns:
                filtered_df['charge'] = filtered_df['scan_number'].map(lambda x: spectra_dict.get(x, {}).get('charge'))    
        
            try: # If MBR was run, the mz is already in the predictions file
                filtered_df['meas_mz'] = filtered_df['mz']
            except KeyError:
                filtered_df['rt'] = filtered_df['scan_number'].map(lambda x: spectra_dict.get(x, {}).get('rt'))
                filtered_df['meas_mz'] = filtered_df['scan_number'].map(lambda x: spectra_dict.get(x, {}).get('meas_mz'))    
            filtered_df['meas_mass'] = filtered_df['charge']*filtered_df['meas_mz'] - (filtered_df['charge']-1)*proton_da
            
            # Initialize number of PSMs to be cummulated later
            filtered_df['no_psms'] = [1 for _ in filtered_df['scan_number']]
            filtered_dict = filtered_df.reset_index().to_dict(orient='index') # Convert to dict for iteration

            iterate_spectras(filtered_dict, noise)

        # Calculate the quantification values for each mzML file
        for _, mzml_id  in tqdm(enumerate(self.mzml_file_dict.keys()), total=len(self.mzml_file_dict), desc='Processing quants for each mzML file'):
            spectra_dict = {}
            exp = MSExperiment()
            MzMLFile().load(self.mzml_file_dict[mzml_id]['file path'], exp)
            max_spectra = exp.getNrSpectra() 
            spectra_dict = {i:{'meas_mz':exp.getSpectrum(i).getPrecursors()[0].getMZ() if exp.getSpectrum(i).getPrecursors() != [] else '',
                               'rt': exp.getSpectrum(i).getRT(),
                               'ms_level': exp.getSpectrum(i).getMSLevel(), 
                               'peaks':exp.getSpectrum(i).get_peaks(), 
                               'charge': exp.getSpectrum(i).getPrecursors()[0].getCharge() if exp.getSpectrum(i).getPrecursors() != [] else ''
                               } for i in range(max_spectra)}
            
            noise = get_noise_for_whole_exp(spectra_dict, noise_boundary=noise_boundary,noise_ms_level=noise_ms_level) 

            # Status file
            self.status_noise.append(f"ID: {mzml_id}, noise: {noise}")
            
            get_quants(spectra_dict, mzml_id, noise)
            self.time_flag(f"Finished running mzml file no. {mzml_id}")
        
        # Return the final dict, which is also saved as a class attribute
        return self.psm_dict

    def make_peptide_table(self, psm_dict = dict()) -> dict:
        """ Pivots the quantification data into a wide format, where the abundance for each peptide 
        in the experiment is moved to a column with a corresponding name. The columns are sorted 
        to combine identical peptides with mass charges within a specified tolerance.

        Parameters
        ----------
        psm_dict : dict
            A dict containing the quantification data for the peptides. This dict should
            include columns such as 'sequence', 'charge', and 'abundance' to be pivoted into the wide format.

        Attributes
        ----------
        peptide_dict : dict
            Dictionary containing peptide information.

        Returns
        -------
        peptide_dict : dict
            A dictionary containing the quantification data for the peptides in a wide format.
            The keys are the peptide sequences and the values are dictionaries with the peptide
            data and corresponding abundance values for each experiment.

        Notes
        -----
        - The function combines peptides with mass charge ratios within tolerance.
        - This function also updates the `self.peptide_dict` attribute of the class with the new data.
        
        """

        # Always check if exception was called, so class will gracefull exit if error was raised
        if self.exception_called:
            return
           
        # Convert to Dataframe for filtering
        psm_df = pd.DataFrame.from_dict(psm_dict, orient='index')
        
        # Group by 'preds' and mass within +- ppm_tolerance
        psm_df['group'] = 0
        sequence_groups = psm_df.groupby(by = ['seq_modifications'])
        # Initalize group number
        group_number = 1
        for sequence, group in sequence_groups:
            group.sort_values(by=['meas_mass'], ascending=True, inplace=True)
            prev_mass = 0
            for i, row in group.iterrows():
                if prev_mass == 0:
                    psm_df.at[i, 'group'] = group_number
                else:
                    # Check if the mass difference is within the tolerance
                    if abs(row['meas_mass'] - prev_mass) <= (self.ppm_tolerance/10**6)*prev_mass:
                        psm_df.at[i, 'group'] = group_number
                    else:
                        group_number += 1
                        psm_df.at[i, 'group'] = group_number

                prev_mass = row['meas_mass']
            
            group_number += 1

        # Combine identical PSMs using the group numbers
        drop_indices = [] 
        for _, (_, group) in tqdm(enumerate(psm_df.groupby('group')), 
                                  total=len(psm_df['group'].unique()), 
                                  desc='Combining identical PSMs'):
            for file_id in group['file_id'].unique():      
                group_file = group[group['file_id'] == file_id]
                best_conf = max(group_file['conf'])
                
                if len(group_file.index.tolist()) > 1:
                    group_line = group_file[group_file['conf'] == best_conf].index[0]
                    group_index_remove = group_file.index.tolist()
                    group_index_remove.remove(group_line)
                    drop_indices.extend(group_index_remove)
                else:
                    group_line = group_file.index[0]
            
                # Filter out rows where 'abundance' is NaN
                abundance_val = [val for val in group_file['abundance'].tolist() if not np.isnan(val)]
                
                # Update DataFrame with combined rows
                psm_df.at[group_line, 'no_psms'] = len(group_file.index)
                psm_df.at[group_line, 'abundance'] = sum(abundance_val)/len(abundance_val) if abundance_val else np.nan
                psm_df.at[group_line, 'conf'] = best_conf
                psm_df.at[group_line, 'meas_mass'] = sum(group_file['meas_mass'])/len(group_file.index)

        # Drop the rows that were combined
        psm_df.drop(index=drop_indices, inplace=True)
        
        # Combine identical PSMs across experiments to peptides
        for _, (_, group) in tqdm(enumerate(psm_df.groupby('sequence')), 
                                  total=len(psm_df['seq_modifications'].unique()), 
                                  desc='Differentiating PSMs with identical sequence but not identical overall'):
            if len(group['group'].unique()) == 1:
                continue
            
            # Assigning groups and adding suffixes to the sequences to ensure error-free pivoting
            group_index = 0
            for group_number in group['group'].unique():
                group_file = group[group['group'] == group_number]

                for i, row in group_file.iterrows():
                    psm_df.at[i, 'sequence'] = f"{row['sequence']}_{group_index}"
                
                group_index += 1
        
        for _, (_, group) in tqdm(enumerate(psm_df.groupby('group')), 
                                  total=len(psm_df['group'].unique()), 
                                  desc='Combining identical PSMs across experiments'):   
            # Calculate the average mass, total PSMs, and best confidence for each peptide
            avg_meas_mass = group['meas_mass'].mean()
            total_psms = group['no_psms'].sum()
            best_conf = group['conf'].max()
            
            # Update the rows for each peptide based on the PSM values
            for i, row in group.iterrows():
                psm_df.at[i, 'meas_mass'] = avg_meas_mass
                psm_df.at[i, 'no_psms'] = total_psms
                psm_df.at[i, 'conf'] = best_conf
            
        # Pivot to wide format using sequence as index 
        wide_df = psm_df.pivot(index=['sequence', 'seq_modifications', 'meas_mass', 'conf', 'no_psms'], columns='file_id', values='abundance').reset_index()
        wide_df.columns = wide_df.columns.to_numpy().tolist()[:5] + [f'abundance_{name}' for name in wide_df.columns.to_numpy().tolist()[5:]]
        
        # Convert to dict and update class attribute
        peptide_dict = wide_df.set_index('sequence').to_dict(orient='index')
        self.peptide_dict = peptide_dict

        return peptide_dict

    def normalizer(self, 
                   peptide_dict = dict(), 
                   type='median') -> dict:
        """ Normalizes the quantification data based on the specified method. The default method is 
        'median', where the column with the highest amount of data points is used as the baseline. 
        Other available methods are 'mean' and 'tic' (Total Ion Current).

        Parameters
        ----------
        peptide_dict : dict
            A dictionary containing the quantification data for the peptides in wide format.

        type : str, optional, default='median'
            The normalization method to use. Options are:

            - 'median' : Normalize using the median of the data.
            - 'mean' : Normalize using the mean of the data.
            - 'tic' : Normalize using the Total Ion Current method.

        Returns
        -------
        normalized_dict : dict
            The normalized peptide dict with additional values for the normalized abundance values. 
            The new values are named with a suffix '_normalized' to indicate the normalization.

        Notes
        -----
        - The normalization method adjusts the abundance values to correct for systematic biases.
        - The 'tic' method sums all intensities in a sample and normalizes each peptide's abundance by the total sum.
        - This method also updates the `self.peptide_dict` attribute with the normalized values.
        
        """

        # Always check if exception was called, so class will gracefull exit if error was raised
        if self.exception_called:   
            return
        
        # Convert to DataFrame for easy full column normalization
        df = pd.DataFrame.from_dict(peptide_dict, orient='index')
        
        # Ensure type is case independent
        type = type.lower()

        # Add variable options to status file
        self.variable_status.append(f"Normalization by {'TIC (Total Ion Current)' if type == 'tic' else type}\n")

        self.time_flag(f"Normalization by {'TIC (Total Ion Current)' if type == 'tic' else type}")
        
        print('Normalizing abundances')
        
        # Calculate the normalized values according to the specified method
        columns, counts = [], []
        for column in df:
            if 'abundance' in column: 
                columns.append(column)
                counts.append(df[column].count())

        if type == 'median' or type == 'mean':
            reference_column = columns[counts.index(max(counts))]
            if type == 'median':
                ref_median = df[reference_column].median()
            elif type == 'mean':
                ref_mean = df[reference_column].mean()

            for column in df:
                if 'abundance' in column: 
                    if type == 'median':
                        col_median = df[column].median()
                        factor = ref_median/col_median
                    elif type == 'mean':
                        col_mean = df[column].mean()
                        factor = ref_mean/col_mean
                    df[f"{column}_normalized"] = df[column]*factor
        
        elif type == 'tic':
            for column in df:
                if 'abundance' in column: 
                    col_sum = df[column].sum()
                    df[f"{column}_normalized"] = df[column]/col_sum
        
        # Convert to dict and update class attribute
        normalized_dict = df.to_dict(orient='index')
        self.peptide_dict = normalized_dict
        
        return normalized_dict

    
    def load_proteome(self, proteome_file = str()):
        """ Loads a proteome file (in FASTA format) into a dictionary, where the key is a tuple 
        of the protein's ID and description, and the value is the protein sequence.

        Parameters
        ----------
        proteome_file : str
            The file path to the proteome file in FASTA format. The file should contain protein 
            sequences, with each entry beginning with a '>' symbol followed by the protein's 
            ID and description, and then the sequence itself on the next lines.

        Attributes
        ----------
        proteome_dict : dict
            Dictionary containing information from the proteome file.
            The keys are tuples of the protein's ID and description, and the values are the protein sequences.

        proteome_description_dict : dict
            Dictionary containing the protein descriptions from the proteome file.
            The keys are tuples of the protein's ID and description, and the values are the protein descriptions.

        Returns
        -------
        None
            The function updates the class attributes `self.proteome_dict` and
            `self.proteome_description_dict` with the loaded protein sequences and their
            corresponding descriptions.

        Notes
        -----
        - The FASTA format should be well-structured, where each protein entry starts with a '>' symbol, followed by an identifier and description (separated by spaces), and the sequence appears on the following lines.
        - If the proteome file contains multiple sequences, each will be parsed and stored in the dictionary.

        """

        # Always check if exception was called, so class will gracefull exit if error was raised
        if self.exception_called:
            return
        
        self.time_flag('Loading proteome file')

        # Use SeqIO to parse the FASTA file 
        proteome_dict_untouched = SeqIO.to_dict(SeqIO.parse(proteome_file, 'fasta'))

        # Write the opening of the file to the status file
        with open(proteome_file, 'r') as f:
            prot_lines = sum(1 for _ in f)
        self.files_opened.add(f"-- Proteome file: {proteome_file}\n---- Line count: {prot_lines}\n")

        # Initialize dictionaries to store protein sequences and descriptions
        self.proteome_dict = {}
        self.proteome_description_dict = {}

        # Iterate through the parsed FASTA file and update the dictionaries
        for k, v in proteome_dict_untouched.items():
            # Use nice format for fasta headers from UniProt (sp for SwissProt, tr for TrEMBL), otherwise use the whole header as key
            k = str(k).split('|')[1] if 'sp|' in k or 'tr|' in str(k) else str(k)
            self.proteome_dict[k] = str(v.seq).replace('I','L')
            self.proteome_description_dict[k] = str(v.description)
        
    def compute_alignments(self, peptide_dict = None, cleavage_length = 4) -> dict:
        """ Runs alignment on the peptide sequences against the supplied FASTA file proteome.
        Updates the peptide table with information about the number of proteins each peptide
        aligns to, the corresponding protein accessions, and the protein matches.
        Initializes the `self.protein_dict` attribute to store protein information when for the
        making of the protein table.

        Parameters
        ----------
        peptide_dict : dict or path to csv/xlsx file, optional
            A dictionary representing the peptide table. If not provided, the function uses 
            the class attribute `self.peptide_dict`. Otherwise, it attempts to load the
            peptide table from the specified file path. 

        cleavage_length : int, default=4
            The length of the cleavage site used for alignment. This parameter is used to output the protein position, 
            with cleavage length being the number of amino acids before and after the peptide sequence in the protein.

        Attributes
        ----------
        protein_dict : dict
            Dictionary containing information about the proteins have peptide matches.

        Returns
        -------
        peptide_dict : dict
            The updated peptide table with additional information about protein matches,
            including the number of protein matches, their accessions, and specific matching proteins.
        
        protein_dict : dict
            Dictionary containing information about the proteins that have peptide matches.
            The keys are tuples of the protein's ID and description, and the values are dictionaries with protein data: accession number, description, peptides aligned, alignment coverage and the length of the protein.

        Notes
        -----
        - The alignment is performed by comparing the peptide sequences against the protein sequences in the provided proteome (FASTA file).
        - For each peptide, the number of protein matches, their accessions, and the specific matching proteins are recorded in the `self.peptide_dict`.
        - The FASTA format should contain the protein sequences in a valid format with '>' headers denoting protein IDs/descriptions followed by the sequence.
        
        """

        # Always check if exception was called, so class will gracefull exit if error was raised
        if self.exception_called:
            return
        
        self.time_flag('Running protein alignment')
        
        # Ensure that the proteome file is loaded for alignment
        if not self.proteome_dict:
            self.__exit__(Break,'Please load a proteome file before running this function')
            return

        # Ensure that peptide dict is either loaded in the class or passed as an argument
        if peptide_dict is None:
            if self.peptide_dict:
                peptide_dict = self.peptide_dict
            else:
                self.__exit__(Break,'Please input or load a peptide table before running this function')
                return

        # Allow for the peptide_dict to be a file path to a csv or xlsx file
        if type(peptide_dict) != dict:
            if type(peptide_dict) == str:
                file_type = peptide_dict.split('.')[-1]
                if file_type == 'csv':
                    peptide_df = pd.read_csv(peptide_dict)
                elif file_type == 'xlsx':
                    peptide_df = pd.read_excel(peptide_dict)
            else:
                self.__exit__(Break,'Please input a dictionary, a csv or an excel file for the peptide table')
                return
            
            try:
                peptide_dict = peptide_df.set_index('sequence').to_dict(orient='index')
            except KeyError:
                print('Peptide column not found in the peptide table')
                print('Please make sure the peptide column is named "sequence"')
            except:
                for seq, group in peptide_df.groupby('sequence'):
                    sequence_id = 0
                    for i, _ in group.iterrows():
                        if len(group) > 1:
                            peptide_df.loc[i, 'sequence'] = seq + '_' + str(sequence_id)
                            sequence_id += 1

                peptide_dict = peptide_df.set_index('sequence').to_dict(orient='index')

        # Update the class attribute with the peptide dictionary
        self.peptide_dict = peptide_dict

        # Initialize the protein dictionary to store protein information
        self.protein_dict = {}

        total_peptides = len(self.peptide_dict)

        # Iterate through the peptide dictionary and align each peptide to the proteome
        for _, (peptide_seq, peptide_data) in tqdm(enumerate(self.peptide_dict.items(), start=1), total=total_peptides, desc='Aligning peptides'):           
            proteins = []
            match_list = []
            cleavage = []
            terminal = []

            peptide_seq_stripped = peptide_seq.split('_')[0]

            # Check all proteins in the proteome dictionary for matches to the peptide sequence
            for accession_id, protein_seq in self.proteome_dict.items():
                protein_match = f"{accession_id} "
                
                # Use regex to find all matches of the peptide sequence
                for match in re.finditer(peptide_seq_stripped, protein_seq):
                    start_pos = match.start()  # Zero-based index
                    end_pos = start_pos + len(peptide_seq_stripped) # Zero-based index
                    # Add the match position to the protein match string
                    protein_match += f"[{start_pos + 1}:{end_pos}]"  # 1-based indexing

                    # Check for teminal positions
                    if start_pos < 2:
                        terminal.append('N-terminal')
                    elif end_pos == len(protein_seq):
                        terminal.append('C-terminal')
                    else:
                        terminal.append('')

                    pre_cleavage = ''
                    post_cleavage = ''

                    # Add the cleavage site information
                    if start_pos < cleavage_length:
                        pre_cleavage = protein_seq[:start_pos]
                        for i in range(0, cleavage_length-len(pre_cleavage)):
                            pre_cleavage = '-' + pre_cleavage
                    else:
                        pre_cleavage = protein_seq[start_pos - cleavage_length:start_pos]
                    
                    if end_pos > len(protein_seq) - cleavage_length:
                        post_cleavage = protein_seq[end_pos:]
                        for i in range(0, cleavage_length-len(post_cleavage)):
                            post_cleavage += '-'
                    else:
                        post_cleavage = protein_seq[end_pos:end_pos + cleavage_length]

                    # Append the peptide environment with the cleavage information
                    cleavage.append(f"[{pre_cleavage}].{peptide_seq_stripped}.[{post_cleavage}]")                
                          
                    if accession_id in self.protein_dict:
                        self.protein_dict[accession_id]['peptide_seqs'].append(peptide_seq)
                        self.protein_dict[accession_id]['coverage'].update(range(int(start_pos+1), int(end_pos + 1 )))
                    else:
                        self.protein_dict[accession_id] = {
                            'accession': accession_id,
                            'description': self.proteome_description_dict[accession_id],
                            'peptide_seqs': [peptide_seq], 
                            'coverage': set(range(int(start_pos+1), int(end_pos + 1))), 
                            'protein_length': len(protein_seq)
                        }

                # Only append if matches were found
                if protein_match != f"{accession_id} ":
                    match_list.append(protein_match)
                    proteins.append(accession_id)

            # Update the peptide dictionary with the protein matches information
            self.peptide_dict[peptide_seq].update({
                'no_protein_matches': len(set(proteins)),
                'protein_matches': ';'.join(proteins),
                'protein_locs': ';'.join(match_list),
                'no_psms': peptide_data['no_psms'],
                'seq_environment': ';'.join(cleavage),
                'terminal': ';'.join(terminal)
            })

            # Calculate the average abundance for each peptide
            sum_ab = 0
            count = 0

            for col in peptide_data.keys():
                if 'abundance' in col: 
                    if np.isnan(self.peptide_dict[peptide_seq][col]):
                        continue
                    sum_ab += self.peptide_dict[peptide_seq][col]
                    count += 1

            if count != 0:
                self.peptide_dict[peptide_seq]['avg_abundance'] = sum_ab/count
            else: 
                self.peptide_dict[peptide_seq]['avg_abundance'] = np.nan

        # Update the protein dictionary with the peptide sequences and coverage information
        for prot in self.protein_dict.keys():
            self.protein_dict[prot]['peptide_seqs'] = ';'.join(sorted(self.protein_dict[prot]['peptide_seqs']))
            self.protein_dict[prot]['coverage'] = len(self.protein_dict[prot]['coverage'])

        self.time_flag('Finished protein alignment')

        # Return the updated peptide and protein class attributes
        return self.peptide_dict, self.protein_dict
        
    def make_protein_table(self, 
                           psm_dict = None, 
                           peptide_dict = None, 
                           protein_dict = None, 
                           quant_method = 'mean', 
                           top_n_peptides = 5) -> dict: 
        """ Groups proteins into clusters based on their shared peptides and properties.
        
        This function analyzes the relationship between peptides and proteins, and groups proteins 
        that share similar peptides. The grouping also takes into account quantification methods 
        and the top N peptides per protein.
        
        Parameters
        ----------
        psm_dict : dict, optional   
            A dictionary containing the PSM data. If not provided, the function uses
            the class attribute `self.psm_dict`.
        
        peptide_dict : dict, optional
            A dictionary containing peptide information. If not provided, the function uses 
            the class attribute `self.peptide_dict`.
            
        protein_dict : dict, optional
            A dictionary containing protein information. If not provided, the function uses 
            the class attribute `self.protein_dict`.
            
        quant_method : str, optional
            The protein quantification method to use for the groups. Options include 'mean' (default),
            'median'.
            
        top_n_peptides : int, optional
            The number peptides to include in the protein quantification. Default is 5.

        Attributes
        ----------	
        grouped_protein_dict : dict
            A dictionary containing the grouped protein information, including the principal protein,
            description, protein group, number of peptides, unique peptides, PSMs, and coverage percentage.
            Also updates the class attributes `self.psm_dict`, `self.peptide_dict`, `self.grouped_protein_dict`
            and `self.protein_dict`.

        psm_columns : list
            A list of the output columns and oclumn order for the PSM table.

        peptide_columns : list
            A list of the output columns and column order for the peptide table.
        
        protein_columns : list
            A list of the output columns and column order for the protein table.
        
        Returns
        -------
        grouped_protein_dict : dict
            A dictionary containing the grouped protein information, including the principal protein,
            description, protein group, number of peptides, unique peptides, PSMs, and coverage percentage.
            Also updates the class attributes `self.psm_dict`, `self.peptide_dict`, `self.grouped_protein_dict`
            and `self.protein_dict`.

        Notes
        -----
        - The protein groups are created based on peptide alignment, with quantification methods applied to summarize protein abundances.
        
        """

        # Always check if exception was called, so class will gracefull exit if error was raised
        if self.exception_called:
            return
        
        if not self.proteome_dict:
            self.__exit__(Break,'Please load a proteome file before running this function')
            return

        self.time_flag('Creating protein table')
        
        if peptide_dict is None:
            if self.peptide_dict:
                peptide_dict = self.peptide_dict
            else:
                self.__exit__(Break,'Please input or load a peptide table before running this function')
                return
        else:
            self.peptide_dict = peptide_dict
        
        if protein_dict is None:
            if self.protein_dict:
                protein_dict = self.protein_dict
            else:
                self.__exit__(Break,'Please input or load a protein table before running this function')
                return
        else:
            self.protein_dict = protein_dict
        
        # Add variable options to status file
        self.variable_status.append(f"Protein quantification method: {quant_method}\n")
        self.variable_status.append(f"Number of peptides per protein (N) for quantification: {top_n_peptides}\n")

        # Processing proteins
        for _, (accession_id, protein_data) in tqdm(enumerate(self.protein_dict.items(), start=1), total=len(self.protein_dict), desc='Processing proteins'):

            # Filter peptides that belong to the current protein
            peptide_list = [self.peptide_dict[seq] for seq in protein_data['peptide_seqs'].split(';')]
            
            # Initialize empty list for the top N peptides
            top_n = []

            # Find top_N most abundant peptides for the current protein
            for seq in peptide_list:
                top_n.sort(reverse=True, key=lambda x: x['avg_abundance'])
                    
                if len(top_n) < top_n_peptides:
                    top_n.append(seq)
                elif seq['avg_abundance'] > top_n[-1]['avg_abundance']:
                    top_n[-1] = seq

            # Update the protein dictionary with peptide information
            self.protein_dict[accession_id]['no_peptides'] = len(peptide_list)
            self.protein_dict[accession_id]['no_unique_peptides_protein'] = sum(1 for seq in peptide_list if seq['no_protein_matches'] == 1)
            self.protein_dict[accession_id]['no_psms'] =  sum(seq['no_psms'] for seq in peptide_list)
            self.protein_dict[accession_id]['coverage_pct'] = round(protein_data['coverage']*100/protein_data['protein_length'])
            
            # Calculate the average abundance for the top N peptides for each protein by the specified method
            for col, abundance in peptide_list[0].items():
                if 'abundance' in col and col != 'avg_abundance':
                    if col not in self.protein_dict[accession_id]:
                        self.protein_dict[accession_id][col] = np.nan
                    
                    counter = 0

                    for seq in top_n:
                        if not np.isnan(abundance):
                            counter += 1
                            if np.isnan(self.protein_dict[accession_id][col]):
                                self.protein_dict[accession_id][col] = abundance
                            else:
                                self.protein_dict[accession_id][col] += abundance

                    if quant_method == 'mean':
                        if counter != 0:
                            self.protein_dict[accession_id][col] = self.protein_dict[accession_id][col]/counter
        
        self.time_flag('Processed proteins, starting grouping')

        # Function for creating the principal protein in the grouped table
        def create_principal_protein(prot = str):
            
            if prot not in self.grouped_protein_dict.keys():
                self.grouped_protein_dict[prot] = {
                    'principal_protein': prot,
                    'description': self.protein_dict[prot]['description'],
                    'protein_group': set(),
                    'no_peptides': 0,
                    'no_unique_peptides_protein': 0,
                    'no_psms': 0,
                    'peptide_seqs': set(self.protein_dict[prot]['peptide_seqs'].split(';')),
                    'coverage_pct': self.protein_dict[prot]['coverage_pct'], 
                    'protein_length': self.protein_dict[prot]['protein_length']
                }

                for col in self.protein_dict[prot].keys():
                    if 'abundance' in col:
                        self.grouped_protein_dict[prot][col] = self.protein_dict[prot][col]
        
        # Initialize the dict for protein inference table
        self.grouped_protein_dict = {}

        # Iterate over peptides to find the principal protein for each peptide
        for _, peptide_values in tqdm(enumerate(self.peptide_dict.values(), start=1), 
                                      total=len(self.peptide_dict), 
                                      desc='Finding principal proteins'):
            
            # Find all proteins that align to the peptide
            proteins = [elm.split(' ')[0] for elm in peptide_values['protein_locs'].split(';')]
            protein_match_dict = {elm.split(' ')[0]: elm for elm in peptide_values['protein_locs'].split(';')}
    
            # Is there only 1 protein?
            if len(proteins) == 1:
                if proteins[0] == '':
                    continue
                peptide_values['principal_protein'] = proteins[0]
                peptide_values['principal_protein_loc'] = protein_match_dict[proteins[0]]
                create_principal_protein(proteins[0])
                continue
                
            # Does any of the proteins have unique peptides? 
            prot_unique = []
            for prot in proteins:
                if self.protein_dict[prot]['no_unique_peptides_protein'] > 0:
                    prot_unique.append(prot)
            
            if len(prot_unique) != 0:
                peptide_values['principal_protein'] = ';'.join(prot_unique)
                peptide_values['principal_protein_loc'] = ';'.join([protein_match_dict[prot] for prot in prot_unique])
                for prot in prot_unique:
                    create_principal_protein(prot)
                continue
            
            # No natural principal protein exists, so find the protein(s) with the most amount of peptides
            best_nr_protein = []
            max_peptides = -1
            for prot in proteins:
                if self.protein_dict[prot]['no_peptides'] > max_peptides:
                    max_peptides = self.protein_dict[prot]['no_peptides']
                    best_nr_protein = [prot]
                elif self.protein_dict[prot]['no_peptides'] == max_peptides:
                    best_nr_protein.append(prot)
            
            # Is there one "best"?
            if len(best_nr_protein) == 1:
                peptide_values['principal_protein'] = best_nr_protein[0]
                create_principal_protein(best_nr_protein[0])
                continue

            # Of the ones with the most peptides, which has highest coverage?
            best_coverage = []
            max_coverage = -1
            for prot in best_nr_protein:
                if self.protein_dict[prot]['coverage_pct'] > max_coverage:
                    max_coverage = self.protein_dict[prot]['coverage_pct']
                    best_coverage = [prot]
                elif self.protein_dict[prot]['coverage_pct'] == max_coverage:
                    best_coverage.append(prot)
            
            # Is there one best?
            if len(best_coverage) == 1:
                peptide_values['principal_protein'] = best_coverage[0]
                peptide_values['principal_protein_loc'] = protein_match_dict[best_coverage[0]]
                create_principal_protein(best_coverage[0])
                continue
            
            # No way to tell them apart? Pick the first one 
            peptide_values['principal_protein'] = best_coverage[0]
            peptide_values['principal_protein_loc'] = protein_match_dict[best_coverage[0]]
            create_principal_protein(best_coverage[0])
        
        self.time_flag('Found principal proteins')

        # Now, all principal protein are created, assign peptides to the principal protein table  
        for _, peptide_values in tqdm(enumerate(self.peptide_dict.values(), start=1), 
                                      total=len(self.peptide_dict), 
                                      desc='Grouping proteins'):
            
            # Dont look at peptides that dont match to a protein
            if peptide_values['protein_matches'].split(';')[0] == '':
                continue

            # For each principal protein that contains the peptide
            for prot1 in peptide_values['principal_protein'].split(';'):
                # If there is more than one protein, compare them
                if len(peptide_values['protein_matches'].split(';')) > 1:
                    for prot2 in peptide_values['protein_matches'].split(';'): 
                        # Don't compare a protein to itself
                        if prot1 == prot2 or prot2 in self.grouped_protein_dict[prot1]['protein_group']:
                            continue
            
                        # Find all of the peptides that the second protein contains                            
                        prot2_peptides = self.protein_dict[prot2]['peptide_seqs'].split(';')
                        
                        # Find the principal proteins that are linked to the second protein
                        principal_proteins = [self.peptide_dict[seq]['principal_protein'].split(';') for seq in prot2_peptides]
                        principal_proteins = set([element for innerList in principal_proteins for element in innerList])

                        # Find the best principal protein to make sure the peptide is assigned to the correct protein group
                        best_principal_protein = []
                        max_peptides_shared = 0

                        for principal_prot in principal_proteins:
                            peptides_shared = len(set(prot2_peptides) & set(self.protein_dict[principal_prot]['peptide_seqs'].split(';')))
                            
                            if max_peptides_shared < peptides_shared:
                                best_principal_protein = [principal_prot]
                                max_peptides_shared = peptides_shared

                            elif max_peptides_shared == peptides_shared:
                                best_principal_protein.append(principal_prot)
                        
                        for best_principal in best_principal_protein:
                            if best_principal != prot2 and prot2 not in self.grouped_protein_dict[best_principal]['protein_group']:
                                self.grouped_protein_dict[best_principal]['protein_group'].add(prot2)
                                for peptide in prot2_peptides:
                                    self.grouped_protein_dict[best_principal]['peptide_seqs'].add(peptide)

                # Add the PSM count to the principal protein from the peptide(s)            
                self.grouped_protein_dict[prot1]['no_psms'] += peptide_values['no_psms']

                # If the peptide is unique to the protein, add to the unique peptide count for that protein
                if len(peptide_values['protein_matches'].split(';')) == 1:   
                    self.grouped_protein_dict[prot1]['no_unique_peptides_protein'] += 1

        # Add the total number of peptides to the principal proteins and the peptides, and the protein group as strings
        protein_group_list = [self.grouped_protein_dict[prot]['protein_group'] for prot in self.grouped_protein_dict.keys()]
        protein_group_list = [group for sublist in protein_group_list for group in sublist]

        # Update the principal proteins with the information about its aligning peptides
        for prot in self.grouped_protein_dict.keys():
            self.grouped_protein_dict[prot]['no_peptides'] = len(self.grouped_protein_dict[prot]['peptide_seqs'])
            self.grouped_protein_dict[prot]['peptide_seqs'] = ';'.join(sorted(self.grouped_protein_dict[prot]['peptide_seqs']))
            self.grouped_protein_dict[prot]['protein_group'] = f"{prot};{';'.join(sorted(self.grouped_protein_dict[prot]['protein_group']))}"
            self.grouped_protein_dict[prot]['no_protein_groups'] = np.count_nonzero(protein_group_list == prot) + 1

        self.time_flag('Grouped proteins')

        # All peptides in the group
        peptide_group_list = [self.grouped_protein_dict[prot]['peptide_seqs'].split(';') for prot in self.grouped_protein_dict.keys()]
        peptide_group_list = [group for sublist in peptide_group_list for group in sublist]

        # Finding locations for the principal protein
        for sequence, values in self.peptide_dict.items():
            self.peptide_dict[sequence]['sequence'] = sequence.split('_')[0]
            self.peptide_dict[sequence]['no_protein_groups'] = sum([1 for peptide in peptide_group_list if peptide == sequence])
            try:
                start_list = [re.sub(r":\w*\]", "]", prot) for prot in values['principal_protein_loc'].split(';')]
                end_list = [re.sub(r"\[\w*:", "[", prot) for prot in values['principal_protein_loc'].split(';')]
                
                self.peptide_dict[sequence]['principal_protein_start'] = ';'.join(start_list)
                self.peptide_dict[sequence]['principal_protein_end'] = ';'.join(end_list)
            except KeyError:
                continue        
        
        # Update the grouped protein dictionary with the number of unique peptides in the group
        for values in self.grouped_protein_dict.values():
            values['no_unique_peptides_group'] = sum(1 for peptide in values['peptide_seqs'].split(';') if self.peptide_dict[peptide]['no_protein_groups'] == 1)
            values['peptide_seqs'] = ';'.join([peptide.split('_')[0] for peptide in values['peptide_seqs'].split(';')])

        # Update PSM table with protein matches
        if psm_dict is not None:
            if type(psm_dict) == dict:
                self.psm_dict = psm_dict
            else: 
                print('Please input a dictionary for the PSM table')
                return
        
        try:
            for index, values in self.psm_dict.items():
                sequence = values['sequence']

                try:
                    self.peptide_dict[sequence]
                except KeyError:
                    # Check for the sequence with the _0 suffix, which is used for duplicate sequences which are not identical
                    sequence += '_0'
                
                try:
                    if self.peptide_dict[sequence]['protein_matches'].split(';')[0] == '':
                        self.psm_dict[index]['no_protein_matches'] = 0
                        self.psm_dict[index]['principal_protein'] = self.empty_values
                        self.psm_dict[index]['protein_matches'] = self.empty_values
                        self.psm_dict[index]['no_protein_groups'] = self.empty_values
                    else:
                        self.psm_dict[index]['no_protein_matches'] = self.peptide_dict[sequence]['no_protein_matches']
                        self.psm_dict[index]['principal_protein'] = self.peptide_dict[sequence]['principal_protein']
                        self.psm_dict[index]['protein_matches'] = self.peptide_dict[sequence]['protein_matches']
                        self.psm_dict[index]['no_protein_groups'] = self.peptide_dict[sequence]['no_protein_groups']
                except KeyError:
                    print(f"Peptide {sequence.replace('_0','')} not found in peptide table")

        except AttributeError:
            print('No PSM data found')

        # Define the desired column order for each DataFrame (for the output files)
        self.psm_columns = ['sequence', 'seq_modifications', 'modifications', 'principal_protein', 'protein_matches', 'no_protein_matches', 'no_protein_groups', 'no_psms', 'conf', 'abundance', 'charge', 'mz', 'calc_mass', 'meas_mass', 'rt', 'id_scan', 'file_id', 'exp_file']
        
        self.peptide_columns = ['sequence', 'seq_environment', 'seq_modifications', 'terminal', 'principal_protein', 'principal_protein_start', 'principal_protein_end', 'protein_matches', 'protein_locs', 'no_protein_matches', 'no_protein_groups', 'no_psms', 'conf', 'meas_mass']

        peptide_col_ab = []
        peptide_col_norm = []

        # Append columns to the peptide table and sort them in the correct order
        for col in next(iter(self.peptide_dict.values())).keys():    
            if col in self.peptide_columns:
                continue
            if 'abundance' in col:
                if col == 'avg_abundance':
                    continue
                elif 'normalized' in col:
                    peptide_col_norm.append(col)
                else:
                    peptide_col_ab.append(col)
            else:
                if col not in self.peptide_columns:
                    print('Column not found in peptide table:')
                    print(col)

        self.peptide_columns += sorted(peptide_col_ab) + sorted(peptide_col_norm)

        self.protein_columns = ['principal_protein', 'description', 'coverage_pct', 'protein_length', 'protein_group', 'no_protein_groups', 'peptide_seqs', 'no_peptides', 'no_unique_peptides_protein', 'no_unique_peptides_group', 'no_psms'] 
        self.protein_columns += sorted(peptide_col_ab) + sorted(peptide_col_norm)

        # Return the protein inference table class attribute
        return self.grouped_protein_dict

    def write_files(self,  
                    write_protein_table = True, 
                    protein_file_name = str(),
                    write_peptide_table = True, 
                    peptide_file_name = str(),
                    write_psm_table = True, 
                    psm_file_name = str(),
                    write_ungrouped_protein_table = False, 
                    ungrouped_protein_file_name = str(),
                    write_individual_quantification_files = False,
                    individual_quant_file_name = str(),
                    overwrite_all = False):
        """ Writes the output files for the experiment. Options include writing individual abundance files, protein table, 
        peptide table, psm table, and grouped protein table. Requires input for the overwrite statement. If none is given, 
        it will prompt for it.
        
        Parameters
        ----------
        write_protein_table : bool, optional
            Whether to write the protein table to a file. Defaults to True.
            
        protein_file_name : str, optional
            Specific file name for the protein table. Defaults to an empty string, in which case the function uses the 
            default name for the protein file.
            
        write_peptide_table : bool, optional
            Whether to write the peptide table to a file. Defaults to True.
            
        peptide_file_name : str, optional
            Specific file name for the peptide table. Defaults to an empty string, in which case the function uses the 
            default name for the peptide file.
            
        write_psm_table : bool, optional
            Whether to write the psm (Peptide Spectrum Match) table to a file. Defaults to True.
            
        psm_file_name : str, optional
            Specific file name for the psm table. Defaults to an empty string, in which case the function uses the 
            default name for the psm file.
            
        write_ungrouped_protein_table : bool, optional
            Whether to write the ungrouped protein table to a file. Defaults to False.
            
        ungrouped_protein_file_name : str, optional
            Specific file name for the ungrouped protein table. Defaults to an empty string, in which case the function uses the 
            default name for the ungrouped protein file.
            
        write_individual_quantification_files : bool, optional
            Whether to write individual quantification files for each peptide or protein. Defaults to False.
            
        individual_quant_file_name : str, optional
            Specific prefix for individual quantification files. Defaults to an empty string, in which case the function uses 
            the default prefix for the quantification files.
            
        overwrite_all : bool, optional
            Whether to overwrite existing files. If set to True, all existing files will be overwritten without prompt. Defaults 
            to False, in which case the function will prompt the user for confirmation before overwriting any files.
            
        Returns
        -------
        None
            This function does not return any value. It writes the output files for the experiment as specified in the 
            parameters.

        Notes
        -----
        - If the `overwrite_all` parameter is set to True, the function will automatically overwrite any existing files without confirmation.
        - The function handles multiple output files, depending on the specific options provided by the user.
        - If any file names are not specified, default names will be used.
        """

        # This function does not check for the graceful exit as it error handles in itself
        # But if an exit is made and this function is called, any files which can be output will be

        def file_name_check(file = str(), file_name = str()):
            while True:
                try:
                    output_file = open(file_name, 'x')
                    output_file.close()
                    break

                except FileExistsError:
                    if overwrite_all:
                        answ = ''
                    else:
                        answ = input(f"The {file} file already exists.\nTried file path:{file_name}.\nOverwrite? [y/n] ")

                    if answ.lower() == 'y' or overwrite_all == True:
                        print(f"You are overwriting the {file} file.")
                        break
                    else:
                        print(f"You chose not to overwrite the {file} file.")
                        answ = input('Do you wish to input a new filename? [y/n] ')
                        if answ.lower() == 'n':
                            file_name = input('Please input a new filename: ')
                            file_name = self.experiment_path + file_name + '.csv'
                        
                        else:
                            self.__exit__(Break, f"You chose to quit the function, because the {file} file already exists.", None)
                            return
                        
            return file_name
        
        if write_protein_table:
            if protein_file_name == '':
                protein_file_name = f"{self.experiment_path}{self.experiment_name}_protein_table.csv"
            else:
                protein_file_name = self.experiment_path + protein_file_name + '.csv'
            protein_file_name = file_name_check('protein table', protein_file_name)
            
            try:
                df = pd.DataFrame.from_dict(self.grouped_protein_dict, orient='index')
                df.loc[:, df.columns.str.contains('abundance', case=False)] = df.loc[:, df.columns.str.contains('abundance', case=False)].round(2)
                try:
                    df.to_csv(protein_file_name, index=False, columns=self.protein_columns, na_rep=self.empty_values)
                except AttributeError:
                    df.to_csv(protein_file_name, index=False, na_rep=self.empty_values)
                
                self.status_file_creation.append(f"-- Protein table: {protein_file_name}\n---- Line count: {len(self.grouped_protein_dict)}\n")
            
            except AttributeError:
                print('No protein table found')
        
        if write_peptide_table:
            # Remove underscores for sequence if these are present in the peptide table (created for easier grouping)
            for sequence in self.peptide_dict:
                self.peptide_dict[sequence]['sequence'] = sequence.split('_')[0]
            
            if peptide_file_name == '':
                peptide_file_name = f"{self.experiment_path}{self.experiment_name}_peptide_table.csv"
            else:
                peptide_file_name = self.experiment_path + peptide_file_name + '.csv'

            peptide_file_name = file_name_check('peptide table', peptide_file_name)
            
            try:
                df = pd.DataFrame.from_dict(self.peptide_dict, orient='index')
                df.loc[:, df.columns.str.contains('abundance', case=False)] = df.loc[:, df.columns.str.contains('abundance', case=False)].round(2)
                try:
                    df.to_csv(peptide_file_name, index=False, columns=self.peptide_columns, na_rep=self.empty_values)
                except AttributeError:
                    df.to_csv(peptide_file_name, index=False, na_rep=self.empty_values)
                
                self.status_file_creation.append(f"-- Peptide table: {peptide_file_name}\n---- Line count: {len(self.peptide_dict)}\n")
            
            except AttributeError:
                print('No peptide table found')
        
        if write_psm_table:
            if psm_file_name == '':
                psm_file_name = f"{self.experiment_path}{self.experiment_name}_psm_table.csv"
            else:
                psm_file_name = self.experiment_path + psm_file_name + '.csv'

            psm_file_name = file_name_check('PSM table', psm_file_name)
            
            try:
                df = pd.DataFrame.from_dict(self.psm_dict, orient='index')
                df.loc[:, df.columns.str.contains('abundance', case=False)] = df.loc[:, df.columns.str.contains('abundance', case=False)].round(2)
                try:
                    df.to_csv(psm_file_name, index=False, columns=self.psm_columns, na_rep=self.empty_values)
                except AttributeError:
                    df.to_csv(psm_file_name, index=False, na_rep=self.empty_values)
                
                self.status_file_creation.append(f"-- PSM table: {psm_file_name}\n---- Line count: {len(self.psm_dict)}\n")
            
            except AttributeError as e:
                print(e)
                print('No PSM table found')
        
        if write_ungrouped_protein_table:
            if ungrouped_protein_file_name == '':
                ungrouped_protein_file_name = f"{self.experiment_path}{self.experiment_name}_protein_table_ungrouped.csv"
            else:
                ungrouped_protein_file_name = self.experiment_path + ungrouped_protein_file_name + '.csv'

            ungrouped_protein_file_name = file_name_check('ungrouped protein table', ungrouped_protein_file_name)

            try:
                df = pd.DataFrame.from_dict(self.protein_dict, orient='index')
                df.loc[:, df.columns.str.contains('abundance', case=False)] = df.loc[:, df.columns.str.contains('abundance', case=False)].round(2)
                df.to_csv(ungrouped_protein_file_name, index=False, na_rep=self.empty_values)
                
                self.status_file_creation.append(f"-- Protein table ungrouped: {ungrouped_protein_file_name}\n---- Line count: {len(self.protein_dict)}\n")
            
            except AttributeError:
                print('No ungrouped protein table found')
        
        if write_individual_quantification_files:
            for mzml_id in self.mzml_file_dict.keys():
                if individual_quant_file_name == '':
                    individual_quant_file_name_id = f"{self.experiment_path}{self.experiment_name}_quantification_{mzml_id}.csv"
                else:
                    individual_quant_file_name_id = f"{self.experiment_path}{individual_quant_file_name}_{mzml_id}.csv"

                individual_quant_file_name_id = file_name_check(f"Individual quant {mzml_id}", individual_quant_file_name_id)

                try:
                    psm_df = pd.DataFrame.from_dict(self.psm_dict, orient='index')
                    psm_df.loc[:, psm_df.columns.str.contains('abundance', case=False)] = psm_df.loc[:, psm_df.columns.str.contains('abundance', case=False)].round(2)
                    psm_df[psm_df['file_id'].str.contains(mzml_id,case=False)].to_csv(f"{individual_quant_file_name_id}",index=False,na_rep=self.empty_values)
                    self.status_file_creation.append(f"-- Quantification file: {individual_quant_file_name_id}\n---- Line count: {psm_df[psm_df['file_id'].str.contains(mzml_id,case=False)].shape[0]}\n")
                    
                except AttributeError:
                    print('No individual quantification tables found')
              
