from pyopenms import *
import pandas as pd
import numpy as np
import time
import sys

class INQuant():
    def __init__(self,predictions_file, mzML_file_list, script_purpose = None):
        """ 
        It is recommended to always initialize the Quantification class in a with-statement. 
        This will provide a status file upon completion of script with an overview of time usage and computations completed. 
        Predictions must be a single file for each experiment, multiple mzML files for each predictions file is fine. 
        Predictions must contain colums ['spectrum_index', 'targets', 'preds', 'sample'] 
        - their order and additional colums and disregarded and do not affect functionality. 
        Please ensure that the 'sample' column in the predictions file corresponds roughly to the names of mzML files. 
        mzML files must be a list, even if there is only 1.
        If you would like to have the script purpose in the status file at the end, please add that as a string.
        """
        self.filename_dict = {}
        self.exp = MSExperiment()
        self.filename_dict['Predictions'] = predictions_file
        self.filename_dict['mzML file'] = {}
        self.filename_dict['Area file'] = {}
        self.filename_dict['This experiment'] = self.filename_dict['Predictions'].split('/')[-1].replace('_labelled_kpreds.csv','')
        self.filename_dict['Experiment path'] = self.filename_dict['Predictions'].replace(self.filename_dict['Predictions'].split('/')[-1],'')
        self.filename_dict['Status file'] = f"{self.filename_dict['Experiment path']}{self.filename_dict['This experiment']}_status.txt"
        self.fileIDs_dict = {}
        self.mzml_paths_dict = {}
        self.sorted_preds_dict = {}
        self.QuantDictionary = {}
        self.QuantDictionaryIndex = 0
        self.mzMLfileID = str()
        self.noise = float()
        self.tolerance_MZ = 0.05
        self.status_file_creation = []
        self.status_noise = []
        self.timer_counts = []
        self.mzML_file_list = mzML_file_list
        self.files_opened = set()
        if script_purpose:
            self.script_purpose = script_purpose
        else:
            self.script_purpose = ''

        print('Filenames were loaded sucessfully.')

    def __enter__(self):
        """
        Enter function starts the timer and ensures compatibility with with-statements. 
        """
        self.start = time.perf_counter()
        return self
    
    def __exit__(self,*args):
        """
        Exit function writes the status file for the experiment. 
        """
        self.stop = time.perf_counter()
        time_elapsed = self.stop - self.start
        if time_elapsed > 180:
            time_string = f'{time_elapsed/60:.3} minute(s)'
        else:
            time_string = f'{time_elapsed:.4} second(s)'
        status = open(self.filename_dict['Status file'],'w')
        
        status.write(f"Status for script run for experiment: {self.filename_dict['This experiment']}\n\n")
        if self.script_purpose != '':
            status.write(f'Script purpose:\n{self.script_purpose}\n\n')

        status.write(f'Time spent on entire script: {time_string}\n')
        if self.timer_counts != []:
            status.write(f'- Time elapsed for time-flags called underway in this script:\n')
            for flag in self.timer_counts:
                status.write(flag)
        status.write(f"\nMass charge tolerance: {self.tolerance_MZ}\n")
        status.write('\nFiles loaded in this script (excluding mzML):\n')
        status.write(f"-- Predictions file: {self.filename_dict['Predictions']}\n")
        with open(self.filename_dict['Predictions'], "r") as f:
            pred_lines = sum(1 for _ in f)
        status.write(f'---- Line count: {pred_lines}\n')
        
        if len(self.files_opened) != 0:
            for line in self.files_opened:
                status.write(line)

        status.write('\nmzML files loaded in running of this script:\n')
        if self.filename_dict['mzML file'] != {}:
            for key in self.filename_dict['mzML file'].keys():
                status.write(f"ID: {key}, filename: {self.filename_dict['mzML file'][key].split('/')[-1]}\n")
        else:
            status.write('None were loaded.\n')
        
        if self.status_noise != []:
            status.write('\nNoise calculations in this script based on MS2 lowest 5% intensities:\n')
            for elm in self.status_noise:
                status.write(f'{elm}\n')
        
        status.write('\nAll mzML files in this experiment:\n')
        for key in self.fileIDs_dict.keys():
            status.write(f'ID: {self.fileIDs_dict[key]}, filename: {key}\n')
        
        if self.status_file_creation != []:
            status.write(f'\nIn this script, the following data files were created:\n')
            for elm in self.status_file_creation:
                status.write(elm)

        status.close()
        print(f"The script has run, details can be found in {self.filename_dict['Status file']}. Remember to rename the status file if you would like to keep it before running next script in this experiment.")

    def time_flag(self, name):
        """
        If a time_flag is set anywhere within the with-statement, 
        the status file will show how long the script took to reach the time_flag.
        """
        time_stop = time.perf_counter()
        time_elapsed = time_stop - self.start
        if time_elapsed > 180:
            time_string = f'{time_elapsed/60:.3} minute(s)\n'
        else:
            time_string = f'{time_elapsed:.4} second(s)\n'

        self.timer_counts.append(f'-- {name}: {time_string}')
    
    def run(self, write_individual_area_files = False, normalize_abundance = True, output_format = 'Wide', overwrite_quant_file = None):
        """
        The run function collects all functionality of the class into one command, 
        to simpify the process of finding quantification values.
        """
        self.load_preds()

        for i, mzML in enumerate(self.fileIDs_dict):
            self.mzMLfile_clean = self.mzML_file_list[i].split('/')[-1].split('.')[0]
            self.load_mzML(mzML,make_area_files=write_individual_area_files)
            self.time_flag(f'Finished running mzML file no {i+1}')
    
        if len(self.mzML_file_list) > 1:
            self.write_quant_file(normalize=normalize_abundance,output_format=output_format, overwrite=overwrite_quant_file)
        else:
            self.write_quant_file(normalize=False, overwrite=overwrite_quant_file, output_format=output_format)
        
    
    def load_preds(self):
        """Loads predictions into class, and defines file ID's for all mzML files in experiment. 
            If only 1 mzML file is given, fileID will always be 'mzML'.
            """

        def find_ids(list_of_strings):
            names = [string.split('/')[-1].split('.')[0] for string in list_of_strings]
            min_name = min([len(name) for name in names])

            unique_list_slice_start = 0
            unique_list_slice_stop = 0

            for i in range(min_name):
                unique_char = set([string[i] for string in names])
                if len(unique_char) != 1:
                    unique_list_slice_start = i
                    break

            for i in range(1, min_name):
                i *= -1
                unique_char = set([string[i] for string in names])
                if len(unique_char) != 1:
                    unique_list_slice_stop = i+1
                    break

            unique_id = [string[unique_list_slice_start:len(string)+unique_list_slice_stop] for string in names]
            
            return unique_id

        preds = pd.read_csv(self.filename_dict['Predictions'], sep = ',')
        sorted_preds_df = preds[preds['targets'].str.replace('I','L') == preds['preds']].copy()
        sorted_preds_df['spectrum_index'] = sorted_preds_df['spectrum_index']-1 # To match with 0-indexing from mzML files
        
        if len(self.mzML_file_list) > 1:
            samples = sorted_preds_df['sample'].drop_duplicates().to_list()

            
            unique_mzml_id = find_ids(self.mzML_file_list)
            unique_preds_id = find_ids(samples)
            
            for (mzml_index, mzml) in enumerate(unique_mzml_id):
                try: 
                    pred_index = unique_preds_id.index(mzml)
                except ValueError:
                    print('IDs found in predictions file do not match mzML-filenames.')
                    print(unique_mzml_id)
                    print(unique_preds_id)
                    sys.exit('The script terminated.')
                
                filtered_sorted_preds_df = sorted_preds_df[sorted_preds_df['sample'] == samples[pred_index]].copy()

                self.sorted_preds_dict[mzml] = filtered_sorted_preds_df

                clean_name = self.mzML_file_list[mzml_index].split('/')[-1].split('.')[0]
                self.fileIDs_dict[clean_name] = mzml
                self.mzml_paths_dict[clean_name] = self.mzML_file_list[mzml_index]
        
        elif len(self.mzML_file_list) == 1:
            mzml_path = self.mzML_file_list[0]
            mzml = 'mzML'
            filtered_sorted_preds_df = sorted_preds_df.copy()
            self.sorted_preds_dict[mzml] = filtered_sorted_preds_df

            clean_name = mzml_path.split('/')[-1].split('.')[0]
            self.fileIDs_dict[clean_name] = mzml
            self.mzml_paths_dict[clean_name] = mzml_path

        else:
            print("Error, mzML filelist empty. Please input a list of mzML files for quantification")
            sys.exit('The script terminated.')
            
    def load_mzML(self,mzML_id, make_area_files=False):
        """
        Loads mzML file data into memory, 
        and calculates noise for the experiment.
        """

        def get_noise_for_whole_exp(spectra = None):
            """
            Calculates the noise for the experiment mzML file, 
            which will be subtracted from the abundance calculation.
            Noise is calculated as the lowest 5% of the intensities peaks 
            in all MS2 scans for the experiment.
            Noise for each experiment can be found in the status file.
            """
            if spectra is None:
                spec_data = [peak for i in range(self.max_spectra) if self.exp.getSpectrum(i).getMSLevel() == 2 for peak in self.exp.getSpectrum(i).get_peaks()[1]]
            else:
                spec_data = [peak for elm in self.spectra_data.values() if elm['MSLevel'] == 2 for peak in elm['Peaks'][1]]

            # Find lowest 5 percent:
            total_length = len(spec_data)

            amount_for_noise = int(total_length*0.05)
            spec_data.sort()

            spectra_minimums = spec_data[0:amount_for_noise]
            noise = np.mean(spectra_minimums)

            return noise
        
        mzMLfile = self.mzml_paths_dict[mzML_id]

        print('Loading mzML file...')
        MzMLFile().load(mzMLfile, self.exp)
        self.max_spectra = self.exp.getNrSpectra() 

        
        self.spectra_data = {i:{'MZ':self.exp.getSpectrum(i).getPrecursors()[0].getMZ() if self.exp.getSpectrum(i).getPrecursors() != [] else '','RT':self.exp.getSpectrum(i).getRT(),'MSLevel': self.exp.getSpectrum(i).getMSLevel(), 'Peaks':self.exp.getSpectrum(i).get_peaks()} for i in range(self.max_spectra)}

        self.noise = get_noise_for_whole_exp(spectra=self.spectra_data)

        print(f'{self.mzMLfile_clean} loaded.')

        self.mzMLfileID = self.fileIDs_dict[mzML_id]
        self.filename_dict['mzML file'][self.mzMLfileID] = mzMLfile.split('/')[-1].split('.')[0]
        
        area_file = f'{self.filename_dict["mzML file"][self.mzMLfileID]}_areas.csv'
        self.filename_dict['Area file'][self.mzMLfileID] = area_file

        self.status_noise.append(f'ID: {self.mzMLfileID}, noise: {self.noise}')

        self.get_areas(make_file=make_area_files)


    def get_areas(self, make_file = False):
        """ 
        Makes the quantification dataframe with quantification values for the given mzML. 
        If it is not possible to calculate a quantification value, it is passed as NaN. 
        If make_file = True, the function will generate induvidual area file for each mzML file loaded. 
        """

        def area_calc(x_axis, y_axis):
            """Vectorized area calculation."""
            dx = np.diff(x_axis)
            y_min = np.minimum(y_axis[:-1], y_axis[1:])
            y_diff = np.abs(np.diff(y_axis))
            square_area = dx * y_min
            triangle_area = 0.5 * dx * y_diff

            return np.sum(square_area + triangle_area)
        
        def boundaries(intensities, times, RT_detection_scan):
            """
            Defines boundaries for area-calculation for a single identification.
            Boundaries are found from the first local minima, starting from eiter end of the retention time window.
            """
            if len(times) <= 1:
                return

            high_boundary = None
            
            if not times[0]<RT_detection_scan<times[-1]:
                return

            for i in range(len(times)):
                if times[i]> RT_detection_scan:
                    low_boundary = i-1
                    high_boundary = i
                    break
            
            if high_boundary is None:
                return

            for i in range(len(intensities) - 1, high_boundary, -1):
                if intensities[i] < intensities[i - 1]:
                    high_boundary = i
                    break

            for i in range(0, low_boundary):
                if intensities[i] < intensities[i + 1]:
                    low_boundary = i
                    break
            
            if not 0<=low_boundary<high_boundary<len(intensities):
                return

            return low_boundary, high_boundary

        
        def iterate_spectras(sorted_predictions_dict = dict):
            """
            Goes through all relevant spectra for the experiment and calculates abundance for each peptide. 
            """
            
            filtered_specs = {key: spec for key, spec in self.spectra_data.items() if spec['MSLevel'] == 1}
            for key, prediction in sorted_predictions_dict.items():
                target_MZ = prediction['MZ']
                RT_first_scan = self.spectra_data[key]['RT']
                tolerance_RT = 1/100 * RT_first_scan
                
                target_time, target_intensities = [], []
                for index, spec in filtered_specs.items():
                    rt = spec['RT']
                    if abs(rt - RT_first_scan) > tolerance_RT:
                        continue  

                    MZ_values = np.array(spec['Peaks'][0])
                    Intensity_values = np.array(spec['Peaks'][1])
                    
                    mz_mask = np.abs(MZ_values - target_MZ) < self.tolerance_MZ
                    intensity_sum = np.sum(Intensity_values[mz_mask])

                    if intensity_sum > 0:
                        target_intensities.append(intensity_sum)
                        target_time.append(rt)
                
                target_time = np.array(target_time)
                target_intensities = np.array(target_intensities)
                
                if target_time.size == 0:
                    continue
                
                try:
                    low, high = boundaries(target_intensities, target_time, RT_first_scan)
                    fill_target_time = target_time[low:high + 1]
                    fill_target_ints = target_intensities[low:high + 1]
                    spec_area = area_calc(fill_target_time / 60, fill_target_ints - self.noise)
                    spec_area = spec_area
                except (TypeError, ValueError):
                    spec_area = np.NaN
                
                self.QuantDictionary[self.QuantDictionaryIndex] = {'Peptide': prediction['targets'],'m/z': prediction['MZ'],'Experiment': self.mzMLfileID,'Abundance (-noise)': spec_area}
                self.QuantDictionaryIndex +=1

        if not self.filename_dict["mzML file"][self.mzMLfileID]:
            print('Please load a mzML file before running this function')
            return
        
        filtered_df = self.sorted_preds_dict[self.mzMLfileID]
        filtered_df['MZ'] = [self.spectra_data[spec]['MZ'] for spec in filtered_df['spectrum_index']]
        
        drop_indices = []
        for seq, group in filtered_df.groupby('targets'):
            group = group.sort_values(by='MZ')
            for i, base in group.iterrows():
                close_indices = group.index[(abs(group['MZ'] - base['MZ']) < self.tolerance_MZ) & (group.index > i)]
                drop_indices.extend(close_indices)

        filtered_df.drop(index=drop_indices, inplace=True)

        filtered_dict = filtered_df.set_index('spectrum_index').to_dict(orient='index')
        
        iterate_spectras(filtered_dict)

        self.QuantDataframe = pd.DataFrame.from_dict(self.QuantDictionary,orient='index')
        if make_file:
            self.QuantDataframe[self.QuantDataframe['Experiment'].str.contains(self.mzMLfileID,case=False)].to_csv(self.filename_dict['Area file'][self.mzMLfileID],index=False)
            self.status_file_creation.append(f'-- Quantification file: {self.filename_dict["Area file"][self.mzMLfileID]}\n---- Line count: {self.QuantDataframe.shape[0]}\n')
    
    def write_quant_file(self, file = None, output_format = 'Long', normalize = False,overwrite = None):
        """
        This function writes the quantification file for the whole experiment.
        Options include normalizing by median (which will output a file with both original and normalized data),
        pivoting either long or wide format (normalized are always wide),
        and inputting a file path. If no file path included, the file will be located where the Predictions file is located. 
        Requires input for the overwrite statement, if none is given it will prompt for it.
        """

        def normalizer(area_df):
            """
            We normalize by median, using the column with the highest amount of datapoints as the baseline. 
            """
            self.time_flag('Starting normalization')
            columns, counts = [], []
            for column in area_df:
                if 'Abundance' in column: 
                    columns.append(column)
                    counts.append(area_df[column].count())

            reference_column = columns[counts.index(max(counts))]
            ref_med = area_df[reference_column].median()

            for column in area_df:
                if 'Abundance' in column: 
                    col_med = area_df[column].median()
                    factor = ref_med/col_med
                    area_df[f'{column}_normalized'] = area_df[column]*factor
            
            return area_df
        
        def make_wider(df = pd.DataFrame):
            """
            This function pivots the quantification data into a wide format, 
            so that the abundance for each peptide in the 
            experiment is moved to a column with a corresponding name. 
            The columns are sorted, to combine identical peptides with 
            mass charges within tolerance of each other.
            """
            data_dict = df.to_dict(orient='list')
            unique_peptides = set(data_dict['Peptide'])
            
            for seq in unique_peptides:
                indices = [i for i, peptide in enumerate(data_dict['Peptide']) if peptide == seq]
                sorted_indices = sorted(indices, key=lambda x: data_dict['m/z'][x])
                
                while sorted_indices:
                    MZ_list = [data_dict['m/z'][i] for i in sorted_indices]
                    
                    if abs(max(MZ_list) - min(MZ_list)) < self.tolerance_MZ:
                        avg = sum(MZ_list) / len(MZ_list)
                        for i in sorted_indices:
                            data_dict['m/z'][i] = avg
                        break
                    else:
                        group_end = 1
                        for i in range(1, len(MZ_list)):
                            if abs(MZ_list[i] - MZ_list[i - 1]) >= self.tolerance_MZ:
                                break
                            group_end += 1
                        
                        avg = sum(MZ_list[:group_end]) / group_end
                        for i in sorted_indices[:group_end]:
                            data_dict['m/z'][i] = avg
                        
                        sorted_indices = sorted_indices[group_end:]
            
            df_updated = pd.DataFrame(data_dict)
            wide = df_updated.pivot(index=['Peptide', 'm/z'], columns='Experiment', values='Abundance (-noise)').reset_index()
            names_before = wide.columns.to_numpy().tolist()
            new_names = names_before[:2] + [f'{name}: Abundance (-noise)' for name in names_before[2:]]
            wide.columns = new_names
            return wide

        if file is None:
            file = self.filename_dict['Experiment path'] + self.filename_dict['This experiment'] + '_quantification.csv'

        try:
            output_file = open(file, 'x')
            output_file.close()
        
        except FileExistsError:
            if overwrite:
                answ = ''
            else:
                answ = input('The area data file exists, overwrite? [y/n] ')
            if answ == 'y' or overwrite == True:
                print('You are overwriting the area data file.')
            else:
                print('You chose not to overwrite the area data file.')
                return 

        if normalize:
            wide_df = make_wider(self.QuantDataframe)
            output_df = normalizer(wide_df)

        elif output_format.lower() == 'wide':
            output_df = make_wider(self.QuantDataframe)

        elif output_format == 'Long':
            output_df = self.QuantDataframe

        output_df.to_csv(file,index=False)
        self.status_file_creation.append(f'-- Quantification file: {file}\n---- Line count: {output_df.shape[0]}\n')
        return output_df