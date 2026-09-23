###############################################################################
"""

This script is superpose all gpcr crystals with waters to a reference MODEL (Global Alignment)
Then, for each water molecule, a local alignmet is done to refine the position
Alignment are done using PyMol

This is the most cpu demanding process
Meant to use 8 CPUS, adapt it for your resources!

It also generates many files that will be removed at the end of HomolWat process
(up to 110MB in ~5000 files, in folder /REC/ )

IMPRTANT: To be run with "pymol -cqr prep_PDBs.py file.pdb"

Internal water molecules from crystals are taken from HomolWatDB
Homolwat 2 of 6
__author__: Eduardo Mayol
"""
###############################################################################

# import MySQLdb as mdb
from pymol import cmd
import os, sys
import glob
from multiprocessing import Pool
from datetime import datetime


def chunklist(pdb_paths, parts):
    """divide la lista en N listas de similar length"""
    avg = len(pdb_paths) / float(parts)
    out_lists = []
    last = 0.0
    while last < len(pdb_paths):
        out_lists.append(pdb_paths[int(last) : int(last + avg)])
        last += avg
    return out_lists


def refine_water(pdbID, wat, f_out_rmsd_wats):
    """Locally align one crystal water onto the model and save its refined pose.

    Returns True when the water was placed.  A PyMOL failure (normally an empty
    ``template`` selection for waters at the edge of the reference structure)
    returns False so the caller can keep going: one bad water must not abort the
    whole reference set.
    """
    watnum = wat[22:27].strip()
    cmd.do(
        "select wat_env, byres "
        + pdbID
        + " within 10 of (resi "
        + watnum
        + " and resn HOH and "
        + pdbID
        + ")"
    )
    cmd.create("wat_out", "wat_env")
    cmd.do("save " + pdbID + "_wat_" + watnum + ".pdb, wat_out")
    cmd.do("load " + pdbID + "_wat_" + watnum + ".pdb")
    cmd.do("select wat2ref, resi " + watnum + " and " + pdbID + "_wat_" + watnum)
    if cmd.count_atoms("wat2ref") == 0:
        # The water did not survive the save/load round trip.
        return False
    cmd.do("select template, byres wat2ref around 10 and protein")
    rmsd_wat = 99.0
    if cmd.count_atoms("template and name CA") > 0:
        try:
            rmsd_wat = cmd.super(
                pdbID + "_wat_" + watnum + " and name CA", "template and name CA"
            )[0]
        except Exception:
            # `super` needs a sequence match; fall back to a plain alignment.
            try:
                rmsd_wat = cmd.align(
                    pdbID + "_wat_" + watnum + " and name CA", "template and name CA"
                )[0]
            except Exception:
                rmsd_wat = 99.0
    # Always emit both the refined water and its RMSD row: stage 3 expects one
    # of each for every water it collected from the reference.
    info_wat = pdbID + " HOH " + watnum + " " + ("%.3f" % rmsd_wat) + "\n"
    f_out_rmsd_wats.write(info_wat)
    cmd.create("watrefined", "wat2ref")
    cmd.do("save water_" + watnum + "_" + pdbID + ".pdb, watrefined")
    return True


def superimpose(listapdbs):
    """ make global and local alignments for each structure in listapdbs"""
    path_req = listapdbs.pop()
    path_recWat = listapdbs.pop()
    path_out = path_req + "/REC/"
    # Save information of Global RMSD
    f_out_rmsd = open(path_req + "info_rmsds.txt", "a+")

    allres = "ALA+ARG+ASP+GLU+THR+TYR+HIS+HSD+HSE+HSP+HID+HIE+HIP+LEU+ILE+VAL+PRO+GLY+TRP+LYS+PHE+ASN+GLN+MET+CYS+SER+ACE+NH2+DI7+DI8"
    wat_ions = "HOH+NA"
    watPDBs = []
    wat_pdbs_file = open(path_req + "Pdbs_resol_L", "r")
    f_rw = wat_pdbs_file.readlines()
    for line in f_rw:
        recep = line.split("_")
        if recep[1].strip() != []:
            pdbs = recep[1].split(" ")
            for pdb in pdbs:
                if pdb != "" and pdb != "\n":
                    watPDBs.append(pdb.strip())
    allowed_uniprots = []
    filog_tree = open(path_req + "filter_rec_list", "r")
    ord_list = filog_tree.readlines()
    for unip in ord_list:
        allowed_uniprots.append(unip.strip())

    # start PyMol process
    cmd.do("cd " + path_req)
    cmd.do("load protein.pdb")

    for pdb in sorted(listapdbs):
        # FILTER CRYSTALS WITH INTERNAL WATERS (Pdbs_resol_L)
        # Reference files are named <pdbid>_<chain>_<receptor>.pdb upstream;
        # the local wrapper caches them as <pdbid>_<receptor>.pdb.  Accept both.
        stem = pdb[len(path_recWat) : -4]
        fields = stem.split("_")
        pdbname = fields[0]
        unipcode = fields[-1]
        file_chain = fields[1] if len(fields) >= 3 else ""
        pdbID = stem
        pdb_chain = pdbID[:6]  # e.g. "5WQC_A", kept for the RMSD log
        if pdbname in watPDBs and unipcode in allowed_uniprots:
            waters = []
            f_in = open(pdb, "r")
            f_r = f_in.readlines()
            for line in f_r:
                if line.startswith("HETATM"):
                    res_n = line[17:21].strip()
                    if res_n == "HOH":  # or res_n == "NA":
                        waters.append(line)
            # Load crystal PDB; do Super AND Align; keep best
            cmd.do("load " + pdb)
            cmd.do("select prot, resn " + allres + " and " + pdbID)
            cmd.do("select close, (br. prot around 5) and " + pdbID)
            cmd.do(
                "select far, not (byres (close or prot) around 5) and not (close or prot) and "
                + pdbID
            )
            cmd.do("remove far")
            rms_data_super = cmd.super(pdbID + "////CA", "protein////CA")
            rms_data_align = cmd.align(pdbID, "protein")
            rms_super = rms_data_super[0]
            rms_align = rms_data_align[0]
            if file_chain:
                chain = file_chain
            else:
                # No chain encoded in the filename: take it from the
                # reference's own protein CA atoms.
                ref_chains = cmd.get_chains(pdbID + " and name CA")
                chain = ref_chains[0] if ref_chains else ""
            if rms_super > rms_align:
                rmsd_val = "%.3f" % rms_align
                info = pdb_chain + " " + str(rmsd_val) + " ALIGN\n"
            else:
                rmsd_val = "%.3f" % rms_super
                info = pdb_chain + " " + str(rmsd_val) + " " + " SUPER\n"
            f_out_rmsd.write(info)
            if rms_align > rms_super:
                rms_data_super = cmd.super(pdbID, "protein")
            cmd.do("select not resn " + allres + "+" + wat_ions + " and " + pdbID)
            cmd.do("cd " + path_out)
            cmd.do("save " + pdbID + ".pdb" + ", " + pdbID)
            f_out_rmsd_wats = open(path_out + "info_rmsds" + pdbID + "_wats.txt", "w")
            # Start Local alignments and save RMSD value
            for wat in waters:
                watnum = wat[22:27].strip()
                chain_wat = wat[21:22].strip()
                if chain_wat != chain:
                    continue
                try:
                    refine_water(pdbID, wat, f_out_rmsd_wats)
                except Exception as exc:
                    print(
                        "WARNING: water %s of %s could not be refined: %s"
                        % (watnum, pdbID, exc)
                    )


## Determinar variables fijas, paths, etc
req_num = sys.argv[3]  # name of specific folder of the process
path = sys.argv[4]  # path where the above folder is located

path_req = path + req_num + "/"
# Directory holding the crystal structures with resolved internal waters.
# Prefer the explicit argument passed by the caller (argv[5]).  The historical
# `path[:-8] + "REC_WAT/"` derivation only resolves when the jobs directory has
# a particular leaf name; for any other layout it silently globs zero
# references and no waters are ever transferred.
path_recWat = path[:-8] + "REC_WAT/"
if len(sys.argv) > 5 and os.path.isdir(sys.argv[5]):
    path_recWat = sys.argv[5].rstrip("/") + "/"


startTime = datetime.now()


# Specify number of CPUs used.  PyMOL's C state is not fork-safe, so the
# worker pool can lose processes part-way through a run; export
# HOMOLWAT_CPUS=1 to run the alignments sequentially when that happens.
num_cpus = max(1, int(os.environ.get("HOMOLWAT_CPUS", "8")))

folder_RW = glob.glob(path_recWat + "*.pdb")

# print len(folder_RW), "total!"

# if there are >10 pdbs in folder_RW, divide list in N (num_cpus)
# path_req is included in list as pool just accept 1 argument
# inside superimpose function it is removed from list with .pop


if len(folder_RW) > 10:
    div_pathsRW = chunklist(folder_RW, num_cpus)
    for partial_list in div_pathsRW:
        partial_list.append(path_recWat)
        partial_list.append(path_req)
    if num_cpus == 1:
        for partial_list in div_pathsRW:
            superimpose(partial_list)
    else:
        # execute in multiprocess
        pool = Pool(processes=num_cpus)
        try:
            pool.map(superimpose, div_pathsRW)
        finally:
            pool.close()
            pool.join()
else:
    folder_RW.append(path_recWat)
    folder_RW.append(path_req)
    superimpose(folder_RW)

endTime = datetime.now()
print("execution time:", endTime - startTime)

# Finish PyMol process
cmd.do("quit")

