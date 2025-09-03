#!/usr/bin/env python
########################################################################################################################
#                                                                                                                       
# NOTICE WITH SOFTWARE                                                                                                  
#                                                                                                                       
# The Cambridge Crystallographic Data Centre (CCDC) provides various scripts to many users for use with CCDC            
# applications. Some scripts may be library scripts, written at some earlier stage in time and distributed to other     
# users. Other scripts may be written de novo or modified library scripts for distribution to a specific client for a   
# specific purpose. Unless otherwise agreed, CCDC reserves the right to store a modified or de novo script and use      
# that script as part of a library available to other users.                                                            
#
# No warranty: regardless of the intent of the parties, CCDC makes no warranty that any script is fit for any           
# particular purpose.                                                                                                   
#                                                                                                                       
# License grant: By accepting the PROTAC Conformer Generator (PCG) script from CCDC, each user accedes to the following
# terms:
#                                                                                                                       
# - The PCG script remains the property of CCDC. Regardless of any changes made by a user, the original source code and
#   script remains the property of CCDC and users agree to make no claim of ownership thereof.
# - Users are granted a license to use the PCG script for any purpose, and to change or modify (edit) the script to suit
#   specific needs.
# - Users may not share the PCG script (unmodified or modified by the user) with any third party without permission from
#   CCDC.
# - Users will acknowledge the original authors when using PCG and derived scripts in their research.
#                                                                                                                       
# Please note, this script is provided as-is, but is not formally supported by CCDC at this time.                       
#                                                                                                                       
# Originally created by Jason Cole (cole@ccdc.cam.ac.uk) and the Cambridge Crystallographic Data Centre.
# Parallel implementation, bug fixes, and code cleanup by Kepa K. Burusco-Goni (kburusco@ccdc.cam.ac.uk).
# Documentation, bug fixes, and code cleanup by Fabio Montisci (fmontisci@ccdc.cam.ac.uk).
#
########################################################################################################################

import argparse
import logging
import pathlib
import random
import os
import sys
import pickle
import pandas as pd
import shutil

from multiprocessing import cpu_count
from multiprocessing import Pool
from multiprocessing import Manager

from ccdc.conformer import ConformerGenerator, DefaultConformerParameterFileLocator
from ccdc.descriptors import MolecularDescriptors
from ccdc.io import MoleculeReader, MoleculeWriter
from ccdc.molecule import Molecule
from ccdc.protein import Protein
from ccdc.search import SubstructureSearch, MoleculeSubstructure


class ProteinScorer(object):
    """ Assigns scores to complexes of the conformer with the monomer proteins overlaid onto the PROTAC. """

    def __init__(self, protein, args):
        self._protein = protein
        self._searcher = MolecularDescriptors.AtomDistanceSearch(protein)
        self.verbose = args.verbose
        # ClashScore parameters
        self.cs_min_d = args.cs_min_d
        self.cs_max_d = args.cs_max_d
        self.cs_max_height = args.cs_max_height
        # SurfaceScore parameters
        self.ss_min_d = args.ss_min_d
        self.ss_max_d = args.ss_max_d
        self.ss_max_height = args.ss_max_height

        # Set up a neighbour lookup for atoms to filter out 1-2 and 1-3 contacts
        self._neighbour_lookup = {}
        for atom in self._protein.atoms:
            if atom.atomic_number != 1:
                atoms2 = list(atom.neighbours)  # 1-2 contacts
                atoms3 = set()
                for next_atom in atoms2:
                    atoms3.update(set(next_atom.neighbours))  # 1-3 contacts
                self._neighbour_lookup[atom] = list(atoms3)

    def get_num_atoms(self):
        return len(self._protein.atoms)

    def interacting_count(self):
        """
        Calculates Surface Score by counting heavy atom-atom contacts in the range [self.ss_min_d,self.ss_max_d)

        :param self.ss_max_d: the longest distance to include. Default = 4.0 Angstroms
        :param self.ss_min_d: the shortest distance to include. Default = 2.0 Angstroms
        :param self.verbose: print out details. Default = no
        :return: number of heavy atom-atom contacts in the range [self.ss_min_d,self.ss_max_d]
        """
        n_contacts = 0
        atoms = self._protein.atoms
        for atom in atoms:
            if atom.atomic_number != 1:

                immediate_neighbours = self._neighbour_lookup[atom]
                v = len([a for a in self._searcher.atoms_within_range(atom.coordinates, self.ss_max_d) if
                         a not in immediate_neighbours]) - \
                    len([a for a in self._searcher.atoms_within_range(atom.coordinates, self.ss_min_d) if
                         a not in immediate_neighbours])

                if self.verbose == 'yes':
                    print(f"Atom {atom.label:>6s} in {atom.residue_label:>8s} forms {v:>6d} contacts")

                n_contacts += v

        return n_contacts

    def heavy_atom_clash_score(self):
        """
        Calculates Clash Score by counting heavy atom-atom contacts in the range [self.cs_min_d,self.cs_max_d)

        :param self.cs_max_d: the longest distance to include. Default = 2.0 Angstroms
        :param self.cs_min_d: the shortest distance to include. Default = 1.0 Angstroms
        :param self.cs_max_height: maximum height to be considered in the clash score. Default = 10.0
        :param self.verbose: print out details. Default = no
        :return: score and number of atom-atom contacts in the range [self.cs_min_d,self.cs_max_d)
        """
        score = 0.0
        n_contacts = 0
        diff = self.cs_max_d - self.cs_min_d
        atoms = self._protein.atoms

        for atom in atoms:
            if atom.atomic_number != 1:
                immediate_neighbours = self._neighbour_lookup[atom]
                atoms_within_distance = self._searcher.atoms_within_range(atom.coordinates, self.cs_max_d)
                for close_atom in atoms_within_distance:
                    if close_atom not in immediate_neighbours and close_atom.atomic_number != 1:
                        n_contacts += 1
                        d = MolecularDescriptors.atom_distance(close_atom, atom)
                        score_term = 0.0
                        if d <= self.cs_min_d:
                            score_term = self.cs_max_height
                        else:
                            score_term = self.cs_max_height * (self.cs_max_d - d) / diff

                        if self.verbose == 'yes':
                            print(
                                f"Atom {atom.label:>6s} in {atom.residue_label:>8s} and atom {close_atom.label:>6s} in "
                                f"{close_atom.residue_label:>8s} score {score_term:>10.4f} at distance {d:>10.4f}")

                        score += score_term
        return score, n_contacts


class ProteinAligner(object):
    '''
    Align proteins onto a PROTAC conformer to build a ternary complex model.
    The ligands from the monomer protein-ligand complexes are overlaid onto their maximum common
    substructures in the PROTAC. The same transformation matrix from the overlay is then applied
    to the respective proteins.
    '''

    def __init__(self, protac_conformation):
        self._protac_conformation = protac_conformation
        self._molecule = Molecule(identifier="full model")
        self._molecule.add_molecule(self._protac_conformation.copy())

    @property
    def full_model(self):
        return Protein(identifier="full model", _molecule=self._molecule._molecule)

    def align_protein(self, protein, ligand_id=None):

        def equivalent_atom(mol, at):
            for atom in mol.atoms:
                if at.coordinates == atom.coordinates:
                    return atom
            return None

        protein = protein.copy()

        if ligand_id is None:
            ligand = protein.ligands[0]
        else:
            for l in protein.ligands:
                if l.identifier == ligand_id:
                    ligand = l
                    break

        # Get common substructures
        mcss = MolecularDescriptors.MaximumCommonSubstructure()
        mcss.settings.ignore_hydrogens = True

        atoms, bonds = mcss.search(self._protac_conformation, ligand)

        first_atoms = list(zip(*atoms))[0]
        second_atoms = list(zip(*atoms))[1]

        first_bonds = [(bond.bond_type, bond.atoms[0], bond.atoms[1]) for bond in list(zip(*bonds))[0]]
        second_bonds = [(bond.bond_type, bond.atoms[0], bond.atoms[1]) for bond in list(zip(*bonds))[1]]

        mol1 = Molecule(identifier="first")
        mol1.add_atoms(first_atoms)
        new_first_bonds = [(pr[0], equivalent_atom(mol1, pr[1]), equivalent_atom(mol1, pr[2])) for pr in first_bonds]

        mol1.add_bonds(new_first_bonds)

        mol2 = Molecule(identifier="second")
        mol2.add_atoms(second_atoms)

        new_second_bonds = [(pr[0], equivalent_atom(mol2, pr[1]), equivalent_atom(mol2, pr[2])) for pr in second_bonds]
        mol2.add_bonds(new_second_bonds)

        overlay = MolecularDescriptors.Overlay(mol1, mol2)

        # Transform the whole protein
        protein.remove_ligand(ligand.identifier)
        protein.transform(overlay.transformation)
        self._molecule.add_molecule(protein)


class ProtacConformerGenerator(ConformerGenerator):
    """
    This is an experimental PROTAC conformer generator that only samples the linker in the PROTAC.
    Knowledge of the "warhead" and "E3 binding group" conformations is required in advance.
    """

    def __init__(self, end_points, settings=None, skip_minimisation=False, nthreads=1,
                 parameter_locator=DefaultConformerParameterFileLocator()):
        super().__init__(settings, skip_minimisation, nthreads)
        self._end_points = [end_point.copy() for end_point in end_points]
        self.generator_settings = settings

    def _get_matches(self, model, fragment, max_hits_per_structure=None):
        """ Sets up a substructure search of a fragment into a molecule. """
        searcher = SubstructureSearch()
        ss = MoleculeSubstructure(fragment)
        searcher.add_substructure(ss)
        return searcher.search(model, max_hits_per_structure=max_hits_per_structure)

    def _match_fragment(self, model, fragment):
        """ Matches the fragments to the full molecule according to the results of a substructure search."""
        # this needs recoding to use driving of the conformation using torsion angles.
        try:
            hits = self._get_matches(model, fragment)
        except Exception as e:
            raise RuntimeError(f"Matching of {fragment.identifier} onto {model.identifier} - {str(e)}")

        # This needs thought - what do we do if we get topological symmetry? Does it matter
        if len(hits) == 0:
            logging.error(f"no match of {fragment.identifier} onto {model.identifier}")
            raise RuntimeError(f"no match of {fragment.identifier} onto {model.identifier}")
        elif len(hits) > 1:

            # Find first match that is unlocked (or last match)
            for i in range(len(hits)):
                molecule_atoms = hits[i].match_atoms()
                free = True
                for bond in model.bonds:
                    if bond.atoms[0] in molecule_atoms and bond.atoms[1] in molecule_atoms:
                        if ConformerGenerator.torsion_locked(bond):
                            free = False
                            break
                if free:
                    break
            if not free:
                logging.warning(
                    f"match of {fragment.identifier} onto {model.identifier} - part of fragment was already locked")
            logging.warning(f"ambiguous match of {fragment.identifier} onto {model.identifier} - using hit {i}")
        else:
            molecule_atoms = hits[0].match_atoms()

        return molecule_atoms

    def _match_fragment_and_lock(self, model, fragment):
        """
        Drives the fragments into the same conformation they have in the full molecule and locks their acyclic torsions.
        Current limitations in the underlying CSD Python API prevent to do the same for ring conformations.
        """
        molecule_atoms = self._match_fragment(model, fragment)
        lookup = {f_atom: m_atom for m_atom, f_atom in zip(molecule_atoms, fragment.atoms)}

        for bond in fragment.bonds:
            if bond.is_rotatable:
                f_atom1 = [atom for atom in bond.atoms[0].neighbours if atom not in bond.atoms][0]
                f_atom2 = bond.atoms[0]
                f_atom3 = bond.atoms[1]
                f_atom4 = [atom for atom in bond.atoms[1].neighbours if atom not in bond.atoms][0]
                t = MolecularDescriptors.atom_torsion_angle(f_atom1, f_atom2, f_atom3, f_atom4)
                model.set_torsion_angle(lookup[f_atom1], lookup[f_atom2], lookup[f_atom3], lookup[f_atom4], t)

        for bond in model.bonds:
            if bond.atoms[0] in molecule_atoms and bond.atoms[1] in molecule_atoms:
                ConformerGenerator.lock_torsion(bond)

        return model

    def generate(self, mols):
        """ Generates the conformers while keeping the locked fragments internal coordinates constrained."""
        fixed_mols = []
        for mol in mols:
            model = mol
            for end_point in self._end_points:
                try:
                    model = self._match_fragment_and_lock(model, end_point)
                except RuntimeError as e:
                    logging.warning(
                        f"Locking error for {model.identifier} - {str(e)} - {end_point.identifier} not locked")

            fixed_mols.append(model)

        if len(fixed_mols) > 0:
            return super().generate(fixed_mols)
        else:
            raise RuntimeError("No suitable molecules provided for generation")

    def find_closest(self, original, conformers):
        """
        For validation. Compares the generated conformers to experimentally observed structures.

        :param original: experimentally observed conformation
        :param conformers: generated conformations
        :return: conformer with lowest RMSD, its RMSD, its index, nested list containing RMSD values for each conformer
        """
        best_rmsd = 1e9
        best_conformer = None
        best_index = -1
        rmsd_score_table = [[] for _ in range(len(conformers))]

        print("\n       Find the closest conformer to experimentally observed conformation\n")
        print('       Conformer             RMSD')
        print('       ---------     ------------')

        for i in range(len(conformers)):
            conformer = conformers[i]
            rmsd_and_transformation = MolecularDescriptors.Overlay(original, conformer)
            print(f'       {i + 1:>9d} {rmsd_and_transformation.rmsd:16.6f}')
            rmsd_score_table[i] = [i + 1, rmsd_and_transformation.rmsd]

            if rmsd_and_transformation.rmsd < best_rmsd:
                best_conformer = conformer
                best_rmsd = rmsd_and_transformation.rmsd
                best_index = i

        print(f"\n       The best RMSD is {best_rmsd:.6f} for conformer {best_index + 1}")
        return best_conformer, best_rmsd, best_index, rmsd_score_table


def _get_command_line_arguments(args):
    """
    Parse input parameters from the command line.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("fullmol", type=pathlib.Path,
                        help="3D model of the full protac in mol2 or sdf format")
    parser.add_argument("endpoints", type=pathlib.Path, nargs='+',
                        help="3D models of the two fragments to fix in mol2 or sdf format")
    parser.add_argument("-of", "--output_file", type=pathlib.Path, default=pathlib.Path("protacs_conformers.mol2"),
                        help="Output file for conformers in mol2 or sdf format")
    parser.add_argument("-rnk", "--output_ranking", type=pathlib.Path, default=pathlib.Path("protacs_full_results.csv"),
                        help="Output file for the ranking and results in CSV format")
    parser.add_argument("-mf", "--minimize_first", action="store_true",
                        help="Whether to minimise the PROTAC model before generating conformers (NOTE : will alter "
                             "the head and tail structure in the PROTAC)")
    parser.add_argument("-mc", "--max_conformations", type=int, default=None,
                        help="The maximum number of conformations to generate")
    parser.add_argument("-nt", "--number_of_threads", type=int, default=1,
                        help="The number of threads to use. Default = 1")
    parser.add_argument("-mut", "--max_unusual_torsions", type=int, default=None,
                        help="The number of unusual torsions to allow in a conformation")
    parser.add_argument("-rd", "--randomize_first", action="store_true",
                        help="Whether to randomize the head and tail positions first (for testing purposes)")
    parser.add_argument("-fct", "--find_closest_to",
                        help="Search for the closest conformer to the file specified", default=None)
    parser.add_argument("-pta", "--proteins_to_align",
                        help="Proteins to align onto the conformers generated in mol2 or pdb format. "
                             "If used in combination with -fct it will output a csv file named "
                             "'rmsd_and_score.csv'", default=[], nargs='*')
    parser.add_argument("-v", "--verbose", type=str, choices=['yes', 'no'], default='no',
                        help="Very verbose output for debugging. Requires -pta. Default = no.")
    parser.add_argument("-csr", "--clashscore_range", help="Save only models with a clash score in the specified "
                                                           "range (takes 2 floats)",
                        type=float, default=[-100_000_000.0, 100_000_000.0], nargs=2)
    parser.add_argument("-ssr", "--surfacescore_range", help="Save only models with a surface score in the "
                                                             "specified range (takes 2 integers)",
                        type=int, default=[-100_000_000, 100_000_000], nargs=2)
    parser.add_argument("-rst", "--restart", type=str, default='no',
                        help="Flag to activate the continuation of failed calculation. Options [yes|no].")
    # Arguments for scoring functions
    # 1) ClashScore
    parser.add_argument("-cs_min_d", "--cs_min_d", type=float, default=1.0,
                        help="Lower bound threshold for ClashScore scoring function. Default = 1.0")
    parser.add_argument("-cs_max_d", "--cs_max_d", type=float, default=2.0,
                        help="Upper bound threshold for ClashScore scoring function. Default = 2.0")
    parser.add_argument("-cs_max_h", "--cs_max_height", type=float, default=10.0,
                        help="Maximum score value for a res-res interaction in ClashScore scoring function. Default = 10.0")
    # 2) SurfaceScore
    parser.add_argument("-ss_min_d", "--ss_min_d", type=float, default=2.0,
                        help="Lower bound threshold for SurfaceScore scoring function. Default = 2.0")
    parser.add_argument("-ss_max_d", "--ss_max_d", type=float, default=4.0,
                        help="Upper bound threshold for SurfaceScore scoring function. Default = 4.0")
    parser.add_argument("-ss_max_h", "--ss_max_height", type=float, default=10.0,
                        help="Maximum score value for a res-res interaction in SurfaceScore scoring function. Default = 10.0")
    return parser.parse_args(args)


def _export_args(args):
    """
    Function used when the restart option is invoked.
    Export command line arguments from current run to a pickled file.
    """
    with open(str(pathlib.Path('args.pkl')), "wb") as f:
        pickle.dump(args, f)


def _import_args():
    """
    Function used when the restart option is invoked.
    Read in the arguments from a previous failed run from a pickled file.
    """
    if os.path.isfile(str(pathlib.Path('args.pkl'))):
        with open(str(pathlib.Path('args.pkl')), 'rb') as f:
            args = pickle.load(f)
    else:
        print("\n       File args.pkl not found in path.")
        print("       This file is required to restart the calculation.")
        print("       Place args.pkl file in folder or submit a new calculation.")
        print("       Terminating execution.\n")
        sys.exit()
    return args


def _find_not_completed_steps(counters):
    """
    Function used when the restart option is invoked.
    This function check which files are not found in the folder from a previous
    run and will assume they need to be calculated again.
    """

    def _check_not_file(counter):
        file = pathlib.Path(f'full_model_{counter}.mol2').resolve()
        if not file.is_file():
            return int(counter)

    indices = [_check_not_file(x) for x in counters]
    indices = [index for index in indices if index is not None]
    return indices


def _randomize(mol):
    """
    Only for testing purposes. Randomly rotate the head and tail fragments.
    """

    def _ijk():
        return random.choice([0, 1])

    mol.translate([random.random() * 10.0, random.random() * 10.0, random.random() * 10.0])
    axis = (0, 0, 0)
    while sum(axis) == 0:
        axis = (_ijk(), _ijk(), _ijk())

    mol.rotate(axis, int(random.random() * 180))


def _find_the_closest_conformer_to_reference_structure(args, ptcg, molecules):
    """
    Find the closest conformer to the reference structure. If no reference structure is explicitly provided as
    argument, the initial input PROTAC molecule model will be used are reference instead.
    """
    if args.find_closest_to is not None:
        with MoleculeReader(str(pathlib.Path(args.find_closest_to).resolve())) as close_mol:
            observed_molecule = close_mol[0]
    else:
        # If no reference structure is provided the input PROTAC model is used for RMSD calculations
        with MoleculeReader(str(args.fullmol.resolve())) as close_mol:
            observed_molecule = close_mol[0]

    best_conformer, best_rmsd, best_index, rmsd_score_table = \
        ptcg.find_closest(observed_molecule, molecules)
    # Export rmsd_score_table for use in case of restart
    with open('rmsd_score_table.pkl', 'wb') as f:
        pickle.dump(rmsd_score_table, f)

    with MoleculeWriter(str(pathlib.Path("closest.mol2"))) as out_closest:
        out_closest.write(best_conformer)
    # Save best conformer information to use in case of restart a calculation
    protac_best_conformer = {"best_index": best_index,
                             "best_rmsd": best_rmsd,
                             "bestrmsd_clashscore": 0.0,
                             "bestrmsd_surfscore": 0.0,
                             "bestrmsd_fullmodel": 'bestrmsd_fullmodel'}
    with open(str(pathlib.Path("protac_best_conformer.pkl")), "wb") as outfile:
        pickle.dump(protac_best_conformer, outfile)


def _get_protacs_conformations(args, num_threads):
    """
    For a fresh new calculation, this block of code will generate the
    ensemble of PROTAC conformers
    """
    # Read in the input molecules
    count = 0
    molecules = []
    endpoint_mols = []
    for endpoint in args.endpoints:
        endpoint_mol = MoleculeReader(str(endpoint.resolve()))[0]
        if endpoint is None:
            logging.error(f"Unable to read a molecule from {str(endpoint.resolve())}")
        count += 1
        # For debugging, set the IDs so if exceptions get thrown we know
        # which part is the problem.
        endpoint_mol.identifier = endpoint_mol.identifier + f" - Fragment {count}"

        if args.randomize_first:
            _randomize(endpoint_mol)

        endpoint_mols.append(endpoint_mol)

    skip_minimisation = not args.minimize_first

    # Actually generate the conformers and write them out
    ptcg = ProtacConformerGenerator(endpoint_mols, skip_minimisation=skip_minimisation, nthreads=num_threads)
    ptcg.settings.max_conformers = args.max_conformations
    ptcg.settings.max_unusual_torsions = args.max_unusual_torsions

    model_mol = MoleculeReader(str(args.fullmol.resolve()))
    if model_mol[0] is None:
        logging.error(f"Unable to read a molecule from {str(args.fullmol.resolve())}")

    outpath = args.output_file.resolve()
    if outpath.is_dir():
        outpath = outpath.joinpath('protacs_conformers.mol2').resolve()
    outname = str(pathlib.Path(outpath).stem)

    # Write the PROTACs conformers to the output file:
    # 1) A mol2 file
    f = str(pathlib.Path(outname + '.mol2'))
    with MoleculeWriter(f) as w:
        hitlist = ptcg.generate(model_mol)
        for dataset in hitlist:
            for conf in dataset:
                molecules.append(conf.molecule)
                w.write(molecules[-1])
    # 2) A pickled binary file
    # f = str(pathlib.Path(outname+'.pkl'))
    # with open(str(pathlib.Path(f)), 'wb') as pkl:
    #    pickle.dump(molecules, pkl)

    # For validation, find the closest conformer to the observed structure
    _find_the_closest_conformer_to_reference_structure(args, ptcg, molecules)
    return molecules


def _load_protacs_from_file(args):
    """
    Function used when the restart option is invoked. This function loads
    the PROTAC conformers from a file in the folder from a previous run.
    """
    outpath = args.output_file.resolve()
    if outpath.is_dir():
        outpath = outpath.joinpath('protacs_conformers.mol2').resolve()
    outname = str(pathlib.Path(outpath).stem)
    pkl = str(pathlib.Path(outname + '.pkl'))
    mol2 = str(pathlib.Path(outname + '.mol2'))
    if os.path.isfile(pkl):
        with open(pkl, 'rb') as f:
            molecules = pickle.load(f)
    elif os.path.isfile(mol2):
        # print(f" File {outname}.pkl not found in path")
        # print(f" Searching for file {outname}.mol2")
        molecules = MoleculeReader(mol2)
    else:
        print(" ERROR: Files {outname}.pkl and {outname}.mol2")
        print("        containing the PROTACs not found in path.")
        print("        These files are required to restart the calculation.")
        print("        Place pkl/mol2 files in folder or submit a new calculation.")
        print("        Terminating execution.")
        sys.exit()

    return molecules


def _score_proteins(args):
    """
    Calculate the score of the monomer proteins to correct the score of the complex
    """
    corrections = {}
    proteins_to_align = [Protein.from_file(f) for f in args.proteins_to_align]
    scorers = [ProteinScorer(p, args) for p in proteins_to_align]
    corrections['clash'] = sum([s.heavy_atom_clash_score()[0] for s in scorers])
    corrections['surface'] = sum([s.interacting_count() for s in scorers])

    with open(str(pathlib.Path('corrections.pkl')), "wb") as f:
        pickle.dump(corrections, f)


def _score_ternary_complex(args, aligner, corrections):
    """
    Calculate the scores for the complex, correcting for the monomers.
    """
    scorer = ProteinScorer(aligner.full_model, args)
    numatoms = scorer.get_num_atoms()
    aligned_clash = scorer.heavy_atom_clash_score()[0] - corrections['clash']
    aligned_surface = scorer.interacting_count() - corrections['surface']
    return numatoms, aligned_clash, aligned_surface


def _save_protac_ternary_complex(args, aligner, outpath, aligned_score, aligned_count):
    """
    Saving the full complex models to mol2 files if they fall in selected range of scores
    """
    comment = None
    if (args.clashscore_range[0] < aligned_score < args.clashscore_range[1]) \
            and (args.surfacescore_range[0] < aligned_count < args.surfacescore_range[1]):
        comment = "Yes"
        with MoleculeWriter(str(outpath)) as w:
            w.write(aligner.full_model)
    else:
        comment = 'No'
    return comment


def _align_and_score_ternary_complexes(args, counter, queue1):
    """
    Function to align and overlay the monomer proteins onto the conformer.
    """
    outpath = args.output_file.resolve()
    if outpath.is_dir():
        outpath = outpath.joinpath(f'full_model_{counter}.mol2').resolve()
    else:
        outpath = outpath.parents[0].joinpath(f'full_model_{counter}.mol2').resolve()

    # Loading PROTACs conformations and extracting current conformation
    molecules = _load_protacs_from_file(args)
    mol = molecules[counter - 1]

    aligner = ProteinAligner(mol)  # mol is the conformer

    # Loading proteins_to_align object from pickle
    proteins_to_align = [Protein.from_file(f) for f in args.proteins_to_align]

    if os.path.isfile(str(pathlib.Path('corrections.pkl'))):
        with open(str(pathlib.Path('corrections.pkl')), 'rb') as f:
            corrections = pickle.load(f)
    else:
        print(" ERROR: File corrections.pkl not found")
        print("        Terminating execution")
        sys.exit()

    # Aligns and overlays the monomer proteins onto the conformer
    for protein in proteins_to_align:
        aligner.align_protein(protein)

    # Score ternary complex
    numatoms, aligned_clash, aligned_surface = _score_ternary_complex(args, aligner, corrections)

    # Get RMSD from conformer generator PROTACs file
    RMSD = 99999.9
    if os.path.isfile(str(pathlib.Path('rmsd_score_table.pkl'))):
        with open(str(pathlib.Path('rmsd_score_table.pkl')), 'rb') as f:
            rmsd_score_table = pickle.load(f)
            RMSD = rmsd_score_table[counter - 1][1]

    # Print scores and other information on screen. Messages may be jammed by
    # multiple instances trying to print at the same time on screen.
    fname = os.path.basename(outpath)

    # Saving the full complex models to mol2 files if they fall in selected range of scores
    comment = _save_protac_ternary_complex(args, aligner, outpath, aligned_clash, aligned_surface)

    # Printing scores and rankings on screen
    scoreline = f'       {counter:9d}   {fname:23s}'
    scoreline += f' {RMSD:16.6f}'
    scoreline += f' {aligned_clash:16.6f}'
    scoreline += f' {aligned_surface:16d}'
    scoreline += f' {comment:>12s}'
    print(scoreline)

    # Printing to the queue to export rankings
    protac_scores = {}
    protac_scores['Conformer'] = rmsd_score_table[counter - 1][0]
    protac_scores['RMSD'] = rmsd_score_table[counter - 1][1]
    protac_scores['NumAtoms'] = numatoms
    protac_scores['ClashScore'] = aligned_clash
    protac_scores['SurfaceScore'] = aligned_surface
    protac_scores['File'] = fname
    protac_scores['Saved'] = comment
    queue1.put(protac_scores)

    # Loading protac_best_conformer object from pickle
    if os.path.isfile(str(pathlib.Path('protac_best_conformer.pkl'))):
        with open(str(pathlib.Path('protac_best_conformer.pkl')), 'rb') as f:
            protac_best_conformer = pickle.load(f)

    # Saving the best_index parameters if this is the best index element
    if counter == protac_best_conformer['best_index'] + 1:
        protac_best_conformer['bestrmsd_clashscore'] = aligned_clash
        protac_best_conformer['bestrmsd_surfscore'] = aligned_surface
        protac_best_conformer['bestrmsd_fullmodel'] = f'full_model_{counter}.mol2'
        with open(str(pathlib.Path('protac_best_conformer.pkl')), "wb") as outfile:
            pickle.dump(protac_best_conformer, outfile)


def _collector(queue):
    """
    Take the information from the queues where the pool workers write their
    outputs and write the content to a single file
    """
    full_scores = []
    while not queue.empty():
        scores = queue.get()
        full_scores.append(scores)
        queue.task_done()
    return full_scores


def _print_best_protacs_values(df):
    """
    Find and export the best protacs ternary complexes according to each
    scoring function
    """
    scoring = {'RMSD': True,
               'ClashScore': True,
               'SurfaceScore': True}

    for key, value in scoring.items():
        df.sort_values(by=[key], ascending=value, inplace=True)
        mol2file = df['File'].iloc[0]
        print(f"       Best complex according to {key:<16s} {mol2file:21s} saved as protac_best_conformer_{key}.mol2")
        oldfile = str(pathlib.Path(mol2file))
        newfile = str(pathlib.Path(f'protac_best_conformer_{key}.mol2'))
        shutil.copy(oldfile, newfile)


def _clean_up():
    """
    Remove temporary files after successful calculation.
    """

    def _remove(file):
        if os.path.isfile(str(pathlib.Path(file))):
            os.remove(str(pathlib.Path(file)))

    workdir = os.getcwd()
    files = [f for f in sorted(os.listdir(workdir))
             if os.path.isfile(os.path.join(workdir, f)) and
             (f[-3:] == 'pkl')]

    _ = [_remove(file) for file in files]


def main():
    """
    Main body of the script that calls the argument parser and runs the core workflow.
    """

    print("\n\n\n")
    print("  -------------------------------------------------")
    print("                 PROTACS CALCULATION               ")
    print("  -------------------------------------------------")
    print("\n\n")
    print("\n => 1) Read arguments from command line")

    args = sys.argv[1:]
    args = _get_command_line_arguments(args)
    if args.number_of_threads < 1:
        num_threads = 1
    else:
        num_threads = args.number_of_threads
    restart = args.restart.lower()
    indices = []

    num_cores = cpu_count()
    print(f"\n       Found {num_cores} cores in this machine")

    if restart == 'no':
        print(f"\n       Running using {num_threads} cores")
        # FOR A FRESH NEW RUN:
        # 1) Export command line arguments in case a restart calculation is required
        _export_args(args)
        print("\n => 2) This is a fresh new calculation: Generate PROTAC ensemble of conformers")
        # 2) Generate PROTAC ensemble of conformers
        molecules = _get_protacs_conformations(args, num_threads)
        indices = [i for i, _ in enumerate(molecules, start=1)]

    else:
        # FOR A RESTART CALCULATION:
        # 1) Load command line arguments from a previous failed run
        args = _import_args()
        if args.number_of_threads < 1:
            num_threads = 1
        else:
            num_threads = args.number_of_threads
        print(f"\n       Running using {num_threads} cores")
        args.restart = 'yes'
        print("\n => 2) This is a restart calculation: Load PROTAC ensemble of conformers")
        # 2) Load PROTAC ensemble of conformers
        molecules = _load_protacs_from_file(args)
        indices = [i for i, _ in enumerate(molecules, start=1)]
        indices = _find_not_completed_steps(indices)

    # Build and score complexes of the conformer with the monomer proteins
    # overlaid onto the PROTAC
    if len(args.proteins_to_align) > 0:
        # Calculate the score of the monomer proteins so we can correct the score of the complex
        print("\n => 3) Calculate score of monomer proteins")
        _score_proteins(args)

        print("\n => 4) Align and score ternary complexes (slow)")
        header1 = '\n       Conformer   Model                  '
        header1 += '             RMSD'
        header1 += '       ClashScore'
        header1 += '     SurfaceScore'
        header1 += '     Inside-range'
        print(header1)
        header2 = '       ---------   -----------------------'
        header2 += '     ------------' * 4
        print(header2)

        ###################################################################
        ## PARALLEL CODE
        ###################################################################
        manager = Manager()
        queue1 = manager.Queue()  # Queue to handle the scoring functions
        args_index_queues = [(args, index, queue1) for index in indices]
        with Pool(processes=num_threads, maxtasksperchild=1) as pool:
            pool.starmap(_align_and_score_ternary_complexes, args_index_queues)
        pool.close()
        pool.join()
        # Collecting information from pool workers
        complete_scores = _collector(queue1)
        ###################################################################

        # Export protacs scores to CSV file
        df1 = pd.DataFrame(complete_scores)
        df1.sort_values(by=['RMSD'], ascending=True, inplace=True)

        f1 = str(pathlib.Path(args.output_ranking).stem)
        if restart == 'no':
            f1 = str(pathlib.Path(f1 + '.csv'))
        else:
            f1 = str(pathlib.Path(f1 + '_restart.csv'))
        df1.to_csv(f1, sep=',', index=False, encoding='utf-8')

        if os.path.isfile(str(pathlib.Path('protac_best_conformer.pkl'))):
            with open(str(pathlib.Path('protac_best_conformer.pkl')), 'rb') as f:
                protac_best_conformer = pickle.load(f)
            print("\n => 5) Print scores and data for best RMSD ternary complex.\n")
            print(f"       The lowest (best) RMSD model is: {protac_best_conformer['bestrmsd_fullmodel']}")
            print(f"       RMSD          = {protac_best_conformer['best_rmsd']:20.8f}")
            print(f"       clash_score   = {protac_best_conformer['bestrmsd_clashscore']:20.8f}")
            print(f"       surface_score = {protac_best_conformer['bestrmsd_surfscore']:20d}")

        # Find and export the best protacs ternary complexes according to each
        # scoring function
        print("\n => 6) Collect the best protacs ternary complexes.\n")
        _print_best_protacs_values(df1)

        print("\n => 7) Clean up temporary files.")
        # Remove all temporary pickle files required for a restart calculation
        # (this last step will only be reached if a calculations succeeds)
        _clean_up()
        print("\n => 8) Normal termination.\n")


if __name__ == '__main__':
    main()