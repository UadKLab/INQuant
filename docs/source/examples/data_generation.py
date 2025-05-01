import random
import numpy as np
import pandas as pd
import os
from pyopenms import *
from typing import List

# --- Settings ---
NUM_FILES = 3
NUM_SPECTRA = 5000
NUM_PREDICTIONS = 200
MS1_PROBABILITY = 0.8  # 80% MS1, 20% MS2
NUM_PROTEINS = 100
OUTPUT_DIR = '/'.join(os.path.abspath(__file__).split('/')[:-1])
FASTA_FILENAME = "example_proteome.fasta"

random.seed(42) 
np.random.seed(42)

# Mass table
AA_MASS = {
    'A': 71.03711, 'R': 156.10111, 'N': 114.04293, 'D': 115.02694,
    'C': 103.00919, 'E': 129.04259, 'Q': 128.05858, 'G': 57.02146,
    'H': 137.05891, 'I': 113.08406, 'L': 113.08406, 'K': 128.09496,
    'M': 131.04049, 'F': 147.06841, 'P': 97.05276, 'S': 87.03203,
    'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841,
}
PROTON = 1.007276
WATER = 18.01056

AMINO_ACIDS = list(AA_MASS.keys())

# --- Functions ---

def generate_random_protein(length: int) -> str:
    return ''.join(random.choices(AMINO_ACIDS, k=length))

def generate_peptide_from_protein(protein: str) -> str:
    start = random.randint(0, len(protein) - 8)
    end = min(len(protein), start + random.randint(8, 20))
    return protein[start:end]

def generate_random_peptide(length: int) -> str:
    return ''.join(random.choices(AMINO_ACIDS, k=length))

def calculate_mass(peptide: str) -> float:
    mass = 0.0
    i = 0
    while i < len(peptide):
        aa = peptide[i]
        if i+5 <= len(peptide) and peptide[i+1:i+5] == '(ox)':
            mass += AA_MASS[aa] + 15.9949
            i += 5
        elif i+7 <= len(peptide) and peptide[i+1:i+7] == '(phos)':
            mass += AA_MASS[aa] + 79.9663
            i += 7
        else:
            try:
                mass += AA_MASS[aa]
            except KeyError:
                print(peptide)
            i += 1
    return mass + WATER

def calculate_mz(mass: float, charge: int) -> float:
    return (mass + charge * PROTON) / charge

def create_spectrum(peaks: List[float], sharpness: float, ms_level: int) -> MSSpectrum:
    spectrum = MSSpectrum()
    spectrum.setMSLevel(ms_level)
    mz_array = np.random.uniform(400, 1600, 350).tolist() + [mz + np.random.normal(0, sharpness) for mz in peaks]
    intensity_array = np.random.uniform(500, 2000, 350).tolist() + np.random.uniform(30000, 60000, len(peaks)).tolist()
    indices = np.argsort(mz_array)
    spectrum.set_peaks((np.array(mz_array)[indices], np.array(intensity_array)[indices]))
    return spectrum

def create_ms2_spectrum(precursor_mz: float, charge: int, sharpness: float) -> MSSpectrum:
    spectrum = MSSpectrum()
    spectrum.setMSLevel(2)
    fragments = sorted(precursor_mz * np.random.uniform(0.3, 0.7, 40) + np.random.normal(0, sharpness, 40))
    intensities = np.random.uniform(3000, 7000, 40)
    spectrum.set_peaks((np.array(fragments), np.array(intensities)))
    precursor = Precursor()
    precursor.setMZ(precursor_mz)
    precursor.setCharge(charge)
    spectrum.setPrecursors([precursor])
    return spectrum

def write_fasta(proteins: List[str], filepath: str):
    with open(filepath, "w") as f:
        for idx, sequence in enumerate(proteins):
            f.write(f">protein_{idx+1}\n")
            for i in range(0, len(sequence), 60):
                f.write(sequence[i:i+60] + "\n")

def add_modifications(peptide: str) -> str:
    mod_peptide = ""
    for aa in peptide:
        if aa == 'M' and random.random() < 0.2:
            mod_peptide += 'M(ox)'
        elif aa in ('S', 'T', 'Y') and random.random() < 0:
            mod_peptide += f'{aa}(phos)'
        else:
            mod_peptide += aa
    return mod_peptide

def simulate_protein_expression(proteins: List[str]) -> List[str]:
    expressed_proteins = []
    for protein in proteins:
        if random.random() < 0.7:  # 70% chance of being expressed
            expressed_proteins.append(protein)
    return expressed_proteins

def add_realistic_noise(spectrum: MSSpectrum, noise_level: float = 0.01):
    mz, intensity = spectrum.get_peaks()
    noise = np.random.normal(0, noise_level, len(mz))
    intensity += noise
    spectrum.set_peaks((mz, intensity))

# --- Main ---

random.seed(42)
os.makedirs(OUTPUT_DIR, exist_ok=True)
prediction_rows = []

# Step 1: Generate shared peptides
shared_peptides = [generate_random_peptide(random.randint(8, 15)) for _ in range(75)]

# Step 2: Generate base proteins
proteins = [list(generate_random_protein(random.randint(300, 600))) for _ in range(NUM_PROTEINS)]

# Step 3: Embed shared peptides into multiple proteins
for peptide in shared_peptides:
    num_insertions = random.randint(2, 6)
    chosen_proteins = random.sample(proteins, num_insertions)
    for protein in chosen_proteins:
        insert_pos = random.randint(0, len(protein) - len(peptide))
        protein[insert_pos:insert_pos+len(peptide)] = list(peptide)

# Convert protein lists back to strings
proteins = [''.join(p) for p in proteins]

# Step 4: Simulate expression and write FASTA
expressed_proteins = simulate_protein_expression(proteins)
fasta_path = os.path.join(OUTPUT_DIR, FASTA_FILENAME)
write_fasta(expressed_proteins, fasta_path)

# Step 5: Generate peptides from expressed proteins
proteome_peptides = list(set(generate_peptide_from_protein(p) for p in expressed_proteins for _ in range(3)))
random.shuffle(proteome_peptides)

peptides = []
seen_peptides = set()
while len(peptides) < NUM_PREDICTIONS:
    peptide = proteome_peptides.pop() if len(peptides) < int(NUM_PREDICTIONS * 0.7) else generate_random_peptide(random.randint(8, 20))
    mod_peptide = add_modifications(peptide)
    if mod_peptide not in seen_peptides:
        peptides.append(mod_peptide)
        seen_peptides.add(mod_peptide)

peptides_info = [(pep, calculate_mz(calculate_mass(pep), charge), charge) for pep in peptides for charge in [1, 2, 3]]

# Step 6: Generate mzML files
for file_idx in range(NUM_FILES):
    exp = MSExperiment()
    filename = f"file_{file_idx+1}.mzML"
    scan_number = 1
    rt = 0.0
    chosen_peptides = random.sample(peptides_info, NUM_PREDICTIONS // NUM_FILES)

    for _ in range(NUM_SPECTRA):
        if random.random() < MS1_PROBABILITY or rt == 0:
            peaks = [mz for (_, mz, _) in chosen_peptides if random.random() < 0.7]
            spec = create_spectrum(peaks, random.uniform(0.003, 0.007), 1)
        else:
            peptide, mz, charge = random.choice(chosen_peptides)
            peak_clarity = random.uniform(0.003, 0.006)
            spec = create_ms2_spectrum(mz, charge, peak_clarity)
            log_prob = np.interp(peak_clarity, [0.003, 0.006], [-0.1, -1.0])
            prediction_rows.append({
                "scan_number": scan_number,
                "preds": peptide,
                "log_probs": round(log_prob, 4),
                "file": filename
            })

        rt += random.uniform(0.01, 0.2)
        spec.setRT(rt)
        spec.setNativeID(f'scan={scan_number}')
        add_realistic_noise(spec)
        exp.addSpectrum(spec)
        scan_number += 1

    mzml_path = os.path.join(OUTPUT_DIR, filename)
    MzMLFile().store(mzml_path, exp)

# Step 7: Save predictions CSV
df = pd.DataFrame(prediction_rows)
csv_path = os.path.join(OUTPUT_DIR, "predictions.csv")
df.to_csv(csv_path, index=False)

print(f"Done! Files written to {OUTPUT_DIR}")
