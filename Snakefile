from pathlib import Path

ROSETTA_BIN = "/mnt/speedy/shared/rosetta.source.release-371/main/source/bin"
ROSETTA_DB  = "/mnt/speedy/shared/rosetta.source.release-371/main/database"
RELAX       = f"{ROSETTA_BIN}/relax.default.linuxgccrelease"

ROOT        = "/mnt/speedy/shared/relax_plus_IA_10r_2025_08_11"
INPUT_DIR   = f"{ROOT}/input_data/renamed_2025_08_15"
OUT_RELAX   = f"{ROOT}/output_data/relax"

binding_energies_2025_08_13 = "binding_energies_2025_08_13.csv" 

pdb_files = list(Path(INPUT_DIR).glob("*.pdb"))

MUTANTS = ["job_001_F27R_mutant", "job_002_T28Y_mutant", "job_003_S30L_mutant", "job_004_Y32D_mutant", "job_005_S33A_mutant", "job_006_S54A_mutant", "job_007_S56V_mutant", "job_008_R98L_mutant", "job_009_G99S_mutant", "job_010_A102H_mutant"]

INTERFACE_ANALYZER = f"{ROSETTA_BIN}/InterfaceAnalyzer.linuxgccrelease"
OUT_IA = f"{ROOT}/output_data/interface_2025_08_13"
ia_outputs = [f"{OUT_IA}/{p.stem}_IA.sc" for p in pdb_files]


# Define output paths for relaxed structures and scores
relaxed_outputs = [f"{OUT_RELAX}/{p.stem}_relaxed.pdb" for p in pdb_files]
score_outputs   = [f"{OUT_RELAX}/{p.stem}.sc"           for p in pdb_files]

rule all:
    input:
        relaxed_outputs,
        score_outputs,
        ia_outputs,
        binding_energies_2025_08_13

rule rename_wt:
    input:
        pdb=lambda wc: str(INPUT_DIR / f"{wc.name}.pdb")
    output:
        pdb=f"INPUT_DIR/renamed_2025_08_15/{{name}}.pdb"
    run:
        import os
        os.makedirs(INPUT_DIR / "renamed", exist_ok=True)
        shell("python scripts/rename_wt_chains.py {input.pdb} {output.pdb}")

rule relax:
    threads: 1
    input:
        pdb=lambda wc: f"{INPUT_DIR}/renamed_2025_08_15/{wc.name}.pdb" if "wt" in wc.name.lower() else f"{INPUT_DIR}/{wc.name}.pdb"

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

rule interface_analyzer:
    threads: 1
    input:
        pdb=f"{OUT_RELAX}/{{name}}_relaxed.pdb"
    output:
        score=f"{OUT_IA}/{{name}}_IA.sc",
        log=f"{OUT_IA}/{{name}}_IA.log"
    shell:
        r"""
        mkdir -p "{OUT_IA}"

        "{INTERFACE_ANALYZER}" \
            -database "{ROSETTA_DB}" \
            -s "{input.pdb}" \
            -fixedchains B C \
            -pack_input \
            -pack_separated \
            -out:file:score_only "{output.score}" \
            > "{output.log}" 2>&1
        """

rule parse_interface_scores:
    input:
        sc_files = [f"{OUT_RELAX}/{mut}.sc" for mut in MUTANTS],
        exp_wt = f"{OUT_RELAX}/6att_wt.sc",
        af_wt  = f"{OUT_RELAX}/6att_wt_af_model.sc"
    output:
        csv = binding_energies_2025_08_13
    run:
        import pandas as pd
        import os

        def read_total_score(file_path):
            """Read total_score from a Rosetta .sc file"""
            with open(file_path) as f:
                for line in f:
                    if line.startswith("SCORE: total_score"):
                        headers = line.strip().split()[1:]
                    elif line.startswith("SCORE:") and not line.startswith("SCORE: total_score"):
                        values = line.strip().split()[1:]
                        score_dict = dict(zip(headers, values))
                        try:
                            return float(score_dict["total_score"])
                        except KeyError:
                            raise KeyError(f"'total_score' not found in {file_path}")
            raise ValueError(f"No SCORE lines found in {file_path}")

        # Read wild-type scores
        exp_wt_score = read_total_score(input.exp_wt)
        af_wt_score  = read_total_score(input.af_wt)

        # Collect results
        results = []
        for mut_sc in input.sc_files:
            if not os.path.exists(mut_sc):
                print(f"[WARNING] Missing {mut_sc}, skipping")
                continue
            mut_name = os.path.basename(mut_sc).replace(".sc", "")
            mut_score = read_total_score(mut_sc)
            ddG_exp = mut_score - exp_wt_score
            ddG_af  = mut_score - af_wt_score
            results.append({
                "mutant": mut_name,
                "score": mut_score,
                "ddG_exp": ddG_exp,
                "ddG_af": ddG_af
            })

        # Save CSV
        df = pd.DataFrame(results)
        df.to_csv(output.csv, index=False)

