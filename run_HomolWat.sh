#!/bin/bash -l
#
# HomolWat driver.
#
# All paths are resolved relative to the CURRENT WORKING DIRECTORY, so run it
# from the directory that holds the HomolWat scripts and their data:
#
#     cd /home/ji.cheng4-umw/simulation_pipeline/HomolWat
#     ./run_HomolWat.sh <model.pdb>
#
# Expected in the cwd:
#     prep_all.py, prep_PDBs_mproc.py, ...     the HomolWat scripts
#     Pdbs_resol_ACT, Pdbs_resol_INACT         reference receptor lists
#     scripts/ref_gpcr.*                       BLAST database
#     REC_WAT/                                 reference crystal waters
#     <model.pdb>                              the model to be watered
#
# Output: <cwd>/validation/<XXXX>_HW/<XXXX>_HW.pdb
#
# Environment overrides:
#     PYMOL_BIN    PyMOL executable (default: first pymol found on PATH)
#     PYTHON_BIN   python interpreter (default: python)

set -euo pipefail

start=`date +%s`

# ── Paths: everything is relative to the current working directory ──────────
path="$(pwd)/"                  # HomolWat root
path_scripts="${path}"          # where the .py scripts + Pdbs_resol_* live
path_wats="${path}REC_WAT/"     # reference crystal waters
path_jobs="${path}validation/"  # job output
path_rec="${path}"              # where the MODEL is; reset below from the argument

python_bin="${PYTHON_BIN:-python}"

# PyMOL: honour $PYMOL_BIN, otherwise take the first pymol on PATH
path_pymol="${PYMOL_BIN:-$(command -v pymol || true)}"

# ── Argument + preflight checks ─────────────────────────────────────────────
if [ $# -lt 1 ]; then
    echo "usage: $0 <model.pdb>" >&2
    exit 2
fi

prot="$1"
if [ ! -f "${prot}" ]; then
    echo "ERROR: model '${prot}' not found (paths are relative to $(pwd))" >&2
    exit 2
fi

# prep_all.py opens (path_rec + pdb_in), so keep the model's own directory in
# path_rec and hand it the bare filename.
path_rec="$(cd "$(dirname "${prot}")" && pwd)/"
prot="$(basename "${prot}")"
pdb_folder="${prot:0:4}_HW"
folder_out="${path_jobs}${pdb_folder}"

missing=()
[ -f "${path_scripts}prep_all.py" ]           || missing+=("${path_scripts}prep_all.py")
[ -f "${path_scripts}Pdbs_resol_INACT" ]      || missing+=("${path_scripts}Pdbs_resol_INACT")
[ -f "${path_scripts}scripts/ref_gpcr.pin" ]  || missing+=("${path_scripts}scripts/ref_gpcr.pin (BLAST database)")
[ -d "${path_wats}" ]                         || missing+=("${path_wats}")
[ -n "${path_pymol}" ]                        || missing+=("pymol -- set PYMOL_BIN or add it to PATH")
if [ ${#missing[@]} -gt 0 ]; then
    echo "ERROR: run this from the HomolWat directory; missing:" >&2
    printf '  - %s\n' "${missing[@]}" >&2
    exit 2
fi

mkdir -p "${path_jobs}"

# Let the steps import EM_functions no matter how they are invoked
export PYTHONPATH="${path_scripts}${PYTHONPATH:+:${PYTHONPATH}}"

# prepare folders and MODEL  1/6
"${python_bin}" "${path_scripts}prep_all.py" "${prot}" "${path}" "${path_jobs}" "${path_rec}" "${pdb_folder}" "${path_scripts}"
echo 'running pymol part...'
# Superpose crystals to model 2/6
"${path_pymol}" -cqr "${path_scripts}prep_PDBs_mproc.py" "${pdb_folder}" "${path_jobs}" "${path_wats}" "${path_scripts}" > pymol_info.txt
# Sort and group waters 3/6
"${python_bin}" "${path_scripts}gene_wats_file_refined.py" "${pdb_folder}" "${path_jobs}" "${path_scripts}"

# Chose one of the script to add (or not internal Na ion)
# Add internal water molecules 4/6
"${python_bin}" "${path_scripts}wat_adder_withNA.py" "${pdb_folder}" "${path_jobs}" "${path_scripts}"
#"${python_bin}" "${path_scripts}wat_adder_noNA.py" "${pdb_folder}" "${path_jobs}" "${path_scripts}"

# merge MODEL and waters, renumber waters and prepare output 5/6
"${python_bin}" "${path_scripts}merge_prot_waters.py" "${folder_out}/HW/" "${pdb_folder}"
# prepare pymol sesion with output 6/6
"${path_pymol}" -cqr "${path_scripts}gene_final_pse.py" "${pdb_folder}" "${path_jobs}" "${pdb_folder}"
end=`date +%s`
runtime=$((end-start))
echo "execution time:" $runtime

# Delete files
rm -f "${folder_out}/REC/"*.pdb 2>/dev/null || true

echo "output is located at "${folder_out} "AS" ${pdb_folder}_HW.pdb
