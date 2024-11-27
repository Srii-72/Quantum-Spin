pip install rdkit pyscf streamlit
import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from pyscf import gto, scf
import matplotlib.pyplot as plt

# Streamlit app to input molecule
st.title("Molecule Spin and Position Prediction")

# User input: SMILES string for the molecule
smiles_input = st.text_input("Enter SMILES string (e.g., 'CCO' for ethanol):")

if smiles_input:
    # RDKit: Generate molecule and 3D positions
    mol = Chem.MolFromSmiles(smiles_input)
    if mol is None:
        st.error("Invalid SMILES string. Please check the input.")
    else:
        # Generate 3D structure
        AllChem.EmbedMolecule(mol)
        AllChem.MMFFOptimizeMolecule(mol)

        # Visualize the molecule's 3D structure
        img = Draw.MolToImage(mol)
        st.image(img, caption="Molecule Structure")

        # PySCF: Perform quantum chemistry to predict spin state
        try:
            mol_gto = gto.Mole()
            # For this demonstration, we use a simple approximation to convert SMILES to atoms
            # In a complete solution, a more sophisticated parsing of SMILES would be needed.
            atom_list = Chem.MolToMolBlock(mol).splitlines()[3:-1]  # Extract atom data (rudimentary)
            atoms = []
            for line in atom_list:
                atom_data = line.split()
                atom = atom_data[0]
                coords = list(map(float, atom_data[1:4]))
                atoms.append(f"{atom} {coords[0]} {coords[1]} {coords[2]}")
            
            mol_gto.atom = '\n'.join(atoms)
            mol_gto.basis = 'sto-3g'  # Basis set
            mol_gto.build()

            # Perform Hartree-Fock calculation (Restricted and Unrestricted)
            mf = scf.RHF(mol_gto)
            mf.kernel()

            # Show results: Energy and spin
            st.write("Total Energy (in Hartree):", mf.e_tot)
            st.write("Spin multiplicity:", mf.spin + 1)

            # Unrestricted Hartree-Fock for spin state (open-shell systems)
            mf_uhf = scf.UHF(mol_gto)
            mf_uhf.kernel()

            st.write("Unrestricted Total Energy:", mf_uhf.e_tot)
            st.write("Unrestricted Spin multiplicity:", mf_uhf.spin + 1)

        except Exception as e:
            st.write("Error in quantum chemistry calculation:", e)

Full Framework Implementation
Step 1: Install Necessary Libraries
bash
Copy code
pip install rdkit pyscf streamlit py3Dmol sqlalchemy
•	py3Dmol: For interactive 3D molecular visualization.
•	SQLAlchemy: For database handling.
•	pyscf: For quantum chemistry calculations.
•	rdkit: For molecular manipulation and optimization.
•	streamlit: For building the web interface.
Step 2: Create the Database (SQLAlchemy)
You can use SQLite for local storage or switch to PostgreSQL for a more production-ready solution.
python
Copy code
# db.py - Database setup using SQLAlchemy

from sqlalchemy import create_engine, Column, Integer, String, Float, Text
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker

Base = declarative_base()

class Molecule(Base):
    __tablename__ = 'molecules'
    id = Column(Integer, primary_key=True)
    name = Column(String, unique=True)
    smiles = Column(String)
    energy = Column(Float)
    spin_multiplicity = Column(Integer)
    optimized_structure = Column(Text)

# SQLite database setup (you can replace it with PostgreSQL for a real project)
DATABASE_URL = "sqlite:///molecules.db"
engine = create_engine(DATABASE_URL, echo=True)

Base.metadata.create_all(engine)
Session = sessionmaker(bind=engine)

def get_session():
    return Session()
Step 3: Define the Core Functions (Molecule Prediction)
This module will include the core functions to process molecules (from SMILES input), perform geometry optimization, and predict spin states.
python
Copy code
# molecule_prediction.py

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from pyscf import gto, scf
from db import get_session, Molecule

def generate_3d_structure(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")
    # Generate 3D structure
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol)
    return mol

def quantum_chemistry_prediction(mol):
    # Convert molecule to PySCF molecule object
    mol_gto = gto.Mole()
    atom_list = Chem.MolToMolBlock(mol).splitlines()[3:-1]  # Extract atom data (rudimentary)
    atoms = []
    for line in atom_list:
        atom_data = line.split()
        atom = atom_data[0]
        coords = list(map(float, atom_data[1:4]))
        atoms.append(f"{atom} {coords[0]} {coords[1]} {coords[2]}")
    
    mol_gto.atom = '\n'.join(atoms)
    mol_gto.basis = 'sto-3g'
    mol_gto.build()

    # Perform Hartree-Fock calculation (Restricted and Unrestricted)
    mf = scf.RHF(mol_gto)
    mf.kernel()

    energy = mf.e_tot
    spin_multiplicity = mf.spin + 1

    # Unrestricted Hartree-Fock for spin state (open-shell systems)
    mf_uhf = scf.UHF(mol_gto)
    mf_uhf.kernel()

    unrestricted_energy = mf_uhf.e_tot
    unrestricted_spin_multiplicity = mf_uhf.spin + 1

    return energy, spin_multiplicity, unrestricted_energy, unrestricted_spin_multiplicity

def save_molecule_to_db(smiles: str, energy: float, spin_multiplicity: int, optimized_structure: str):
    session = get_session()
    molecule = Molecule(
        name=smiles,
        smiles=smiles,
        energy=energy,
        spin_multiplicity=spin_multiplicity,
        optimized_structure=optimized_structure
    )
    session.add(molecule)
    session.commit()
    session.close()
Step 4: Streamlit UI for User Input
Here we’ll create the Streamlit application where users can input the molecule and interact with the system.
python
Copy code
# app.py

import streamlit as st
from rdkit.Chem import Draw
from molecule_prediction import generate_3d_structure, quantum_chemistry_prediction, save_molecule_to_db
import py3Dmol

st.title("Molecule Spin and Position Prediction")

# User input: SMILES string for the molecule
smiles_input = st.text_input("Enter SMILES string (e.g., 'CCO' for ethanol):")

if smiles_input:
    try:
        # Step 1: Generate 3D Structure using RDKit
        mol = generate_3d_structure(smiles_input)
        img = Draw.MolToImage(mol)
        st.image(img, caption="Optimized Molecule Structure")

        # Step 2: Perform Quantum Chemistry to Predict Spin and Energy
        energy, spin_multiplicity, unrestricted_energy, unrestricted_spin_multiplicity = quantum_chemistry_prediction(mol)

        st.write(f"Total Energy (in Hartree): {energy}")
        st.write(f"Spin Multiplicity (Restricted): {spin_multiplicity}")
        st.write(f"Unrestricted Energy: {unrestricted_energy}")
        st.write(f"Unrestricted Spin Multiplicity: {unrestricted_spin_multiplicity}")

        # Save the results to the database
        save_molecule_to_db(smiles_input, energy, spin_multiplicity, "Optimized structure here")

        # Step 3: Visualize 3D molecule with py3Dmol
        viewer = py3Dmol.view(width=800, height=400)
        block = Chem.MolToMolBlock(mol)
        viewer.addModel(block, "mol")
        viewer.setStyle({'stick': {}})
        viewer.zoomTo()
        viewer.show()

    except Exception as e:
        st.error(f"Error: {e}")
