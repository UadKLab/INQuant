=========================
Output Data
=========================

Calibration Data
==========================

The ``Calibration`` class outputs calibrated mzML files, which are based on the original mzML files, but with the m/z values corrected. The output files are named according to the original mzML file names, with the suffix "_calibrated" added. For example, if the original mzML file is named "sample.mzML", the calibrated file will be named "sample_calibrated.mzML". The calibrated mzML files are saved in the same directory as the original mzML files unless otherwise specified.
This will be deprecated when the calibration is implemented in the InstaNovo pipeline.

Quantification Data
==========================

The INQuant class can output a variety of data files, based on user specifications. 
Default output files are the PSM table, the peptide table and if protein alignment is computed, the protein table.
Other files, which can be written are the individual quantification files for each mzML file and the ungrouped protein table, which does not include protein inference.

Columns in the output files are as follows:

PSM table
----------------------------

.. list-table:: 
   :header-rows: 1
   :widths: 20 80

   * - Column Name
     - Description
   * - sequence
     - Peptide sequence identified
   * - seq_modifications
     - Sequence including the modifications on specific residues in the sequence
   * - modifications
     - List of modifications present in the peptide
   * - principal_protein*
     - The principal protein which the peptide aligns to
   * - protein_matches*
     - Number of protein entries matching this peptide
   * - no_protein_matches*
     - Number of unique matches to proteins
   * - no_protein_groups*
     - Number of protein groups the peptide is part of
   * - no_psms
     - Number of peptide-spectrum matches (PSMs) combined to this PSM
   * - conf
     - Confidence score or level of identification
   * - abundance
     - Quantitative measure of peptide abundance
   * - charge
     - Charge state of the peptide ion
   * - mz
     - Mass-to-charge ratio of the peptide ion
   * - calc_mass
     - Theoretical calculated mass of the peptide
   * - meas_mass
     - Measured mass from mass spectrometry
   * - rt
     - Retention time during chromatography
   * - id_scan
     - Identification scan number where the peptide was identified by InstaNovo
   * - file_id
     - Identifier for the mzML file containing the scan
   * - exp_file
     - mzML file name

\*only if protein alignment is computed

Peptide table
----------------------------

.. list-table:: 
   :header-rows: 1
   :widths: 20 80

   * - Column Name
     - Description
   * - sequence
     - Peptide sequence identified
   * - seq_environment*
     - Environment of the peptide in the protein sequence
   * - seq_modifications
     - Sequence including the modifications on specific residues in the sequence
   * - terminal*
     - Terminal location, if any
   * - principal_protein*
     - The principal protein which the peptide aligns to
   * - principal_protein_start*
     - Start position of the peptide in the principal protein
   * - principal_protein_end*
     - End position of the peptide in the principal protein
   * - protein_matches*
     - Number of protein entries matching this peptide
   * - protein_locs*
     - Locations of the peptide in the matched proteins
   * - no_protein_matches*
     - Number of unique matches to proteins
   * - no_protein_groups*
     - Number of protein groups the peptide is part of
   * - no_psms
     - Number of peptide-spectrum matches (PSMs) combined to this peptide
   * - conf
     - Confidence score or level of identification
   * - meas_mass
     - Measured mass from mass spectrometry
   * - abundance_[experiment_id]
     - Quantificaton values for each experiment (multiple columns)
   * - abundance_[experiment_id]_normalized
     - Normalized quantification values for each experiment (multiple columns) based on the method specified by the user

\*only if protein alignment is computed

Protein table
----------------------------

.. list-table:: 
   :header-rows: 1
   :widths: 20 80

   * - Column Name
     - Description
   * - principal_protein
     - The principal protein which the peptide aligns to
   * - description
     - Description of the protein
   * - coverage_pct
     - Percentage of the protein sequence covered by identified peptides
   * - protein_length
     - Length of the protein sequence
   * - protein_group
     - Proteins in the group
   * - no_protein_groups
     - Number of protein groups the protein is part of
   * - peptide_seqs
     - List of peptide sequences aligned to this protein
   * - no_peptides
     - Number of peptides aligned to this protein
   * - no_unique_peptides_protein
     - Number of unique peptides aligned to this protein, not shared with other proteins in the group
   * - no_unique_peptides_group
     - Number of unique peptides aligned to this protein group
   * - no_psms
     - Number of peptide-spectrum matches (PSMs) combined to this PSM
   * - abundance_[experiment_id]
     - Quantificaton values for each experiment (multiple columns) for the top N peptides aligning to this protein (N specified by the user)
   * - abundance_[experiment_id]_normalized
     - Normalized quantification values for each experiment (multiple columns) based on the method specified by the user for the top N peptides aligning to this protein (N specified by the user)