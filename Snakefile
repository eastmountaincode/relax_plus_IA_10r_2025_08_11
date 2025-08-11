from pathlib import Path

ROSETTA_BIN = "/mnt/speedy/shared/rosetta.source.release-371/main/source/bin"
ROSETTA_DB  = "/mnt/speedy/shared/rosetta.source.release-371/main/database"
RELAX       = f"{ROSETTA_BIN}/relax.default.linuxgccrelease"

ROOT        = "/mnt/speedy/shared/relax_plus_IA_10r_2025_08_11"
INPUT_DIR   = f"{ROOT}/input_data/pdb"
OUT_RELAX   = f"{ROOT}/output_data/relax"

# Discover all input PDB files
pdb_files = list(Path(INPUT_DIR).glob("*.pdb"))

# Define output paths for relaxed structures and scores
relaxed_outputs = [f"{OUT_RELAX}/{p.stem}_relaxed.pdb" for p in pdb_files]
score_outputs   = [f"{OUT_RELAX}/{p.stem}.sc"           for p in pdb_files]

rule all:
    input:
        relaxed_outputs,
        score_outputs

rule relax:
    threads: 1
    input:
        pdb=lambda wc: f"{INPUT_DIR}/{wc.name}.pdb"
    output:
        pdb=f"{OUT_RELAX}/{{name}}_relaxed.pdb",
        score=f"{OUT_RELAX}/{{name}}.sc",
        log=f"{OUT_RELAX}/{{name}}.log"
    shell:
        r"""
        mkdir -p "{OUT_RELAX}"

        "{RELAX}" \
            -database "{ROSETTA_DB}" \
            -in:file:s "{input.pdb}" \
            -relax:constrain_relax_to_start_coords \
            -nstruct 1 \
            -out::pdb \
            -out:path:pdb "{OUT_RELAX}" \
            -out:suffix _relaxed \
            -out:file:scorefile "{output.score}" \
            > "{output.log}" 2>&1

        mv "{OUT_RELAX}/{wildcards.name}_relaxed_0001.pdb" "{output.pdb}"
        """


