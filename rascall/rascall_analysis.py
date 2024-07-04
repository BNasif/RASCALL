import numpy as np
import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sqlite3
import csv
# Adjust sys.path to include the directory containing plot_NIST.py
from .plot_NIST import NIST_Smile_List, load_NIST_spectra
from rascall.analysis import get_functionals, get_molecules
def calculate_average_properties(functional_dict, functional_code): ## functional_dict is output from get_functionals(), functional code is fg smiles code
    if functional_code not in functional_dict:
        return f"No functional found with code {functional_code}"

    functional = functional_dict[functional_code]
    avg_symmetries = functional.averageSymmetries()

    avg_frequencies = [prop.frequency_average() for symmetry in avg_symmetries for prop in symmetry.properties]
    avg_intensities = [prop.intensity for symmetry in avg_symmetries for prop in symmetry.properties]

    if avg_frequencies and avg_intensities:
        avg_frequency = np.mean(avg_frequencies)
        avg_intensity = np.mean(avg_intensities)
        return avg_frequencies, avg_intensities
    else:
        return "No valid properties found to calculate averages."

def get_list_of_molecule_with_fg(molecule_dictionary, functional, NIST_Smiles):
    molecules_with_test_functional_in_NIST = []
    for molecule_code, molecule_functionals in molecule_dictionary.items():
        if any(functional in s for s in molecule_dictionary.get(molecule_code)):
            if molecule_code in NIST_Smiles:
                molecules_with_test_functional_in_NIST.append(molecule_code)
    if len(molecules_with_test_functional_in_NIST) == 0:
        print ('Requested functional group', functional, 'not in database')
        return

    print ('Number of molecules_with_test_functional in NIST:', len(molecules_with_test_functional_in_NIST))
    return molecules_with_test_functional_in_NIST

def find_and_store_NIST_spectra(molecules_with_functional, functional_dict, functional_to_test, output_csv='nist_spectra_fg.csv'):

    with open(output_csv, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Molecule', 'Wavenumber', 'Absorbance', 'Average Frequency', 'Average Intensity'])  # CSV header
        for molecule in molecules_with_functional:
            try:
                nu, absorb = load_NIST_spectra(molecule, ["wn", "A"], is_smile=True)
                avg_frequencies, avg_intensities, avg_frequency, avg_intensity = calculate_average_properties(functional_dict, functional_to_test)
                for w, a in zip(nu, absorb):
                    writer.writerow([molecule, w, a, avg_frequency, avg_intensity])
            except Exception as e:
                print(f"Error processing molecule {molecule}: {e}")

def main():#Get list of molecules for given functional group. Get band center for fg peak. Find peaks around RASCALL prediction
    NIST_data = NIST_Smile_List()
    NIST_Smiles = NIST_data[0]

    molecule_dictionary = get_molecules()
    functional_dictionary = get_functionals()
    functional = 'CBr'

    molecules_with_functional = get_list_of_molecule_with_fg(molecule_dictionary, functional, NIST_Smiles)
    if molecules_with_functional:
        find_and_store_NIST_spectra(molecules_with_functional, functional, functional_dictionary)


if __name__ == "__main__":
    main()