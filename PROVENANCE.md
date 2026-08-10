# Provenance: `results/all_spliceai_scores_v3.1.tsv`

This document records how the SpliceAI-across-MSA score table was generated, for
publication. It was reconstructed from git history on 2026-08-07; every claim
below is backed by an artefact in this repository and the check that verifies it.

## The data file

| Property | Value |
|---|---|
| Path | `results/all_spliceai_scores_v3.1.tsv` |
| Git blob | `7b372157001460beaf6b0d7db85275e0d98aa27f` |
| First committed | `b4d77bc` "v3.1 analysis", **2026-02-06**, as `results/splice_out_v3.1/all_spliceai_scores_v3.1.tsv` |
| Rows | 4887 (header + 4886 decoy loci) |
| Columns | `decoyID` + 100 species IDs (UCSC 100-way alignment names) |
| Cell format | `[max1][pos1][max2][pos2][d]` |

The working-tree copy is **byte-identical** to the blob first committed on
2026-02-06 (`git hash-object` matches), so the file has not been regenerated or
edited since.

### Cell field meanings

| Field | Meaning |
|---|---|
| `max1` | Max SpliceAI donor score in the canonical 5′ splice-site search window |
| `pos1` | Position (nt from intron start) of `max1` |
| `max2` | Max SpliceAI donor score downstream of `pos1 + 20` — the cryptic/decoy element |
| `pos2` | Position of `max2` — plotted as "loci distance from canonical 5′ splice site" |
| `d` | Normalized intron-length divergence from hg38: `d = abs(n - v) / v`, where `v` is the hg38 ungapped intron length and `n` is this species'. **Not** a phylogenetic distance — it is a length ratio, and it is what drives the length-adjusted `max1` search window (see `find_top2_scores_length_adjusted` in `process_maf_spliceai_v2.5.py`) |

`-1` in any field means the score could not be computed (species absent from the
alignment block, empty/invalid MAF, or all-gap sequence).

## Pipeline that produced the file

```
DecoysInConservedGenes50thpercentile.bed
  │
  ├─ scripts/split_chromosomes.sh              split BED into ~50-locus parts -> chrparts.txt (109 parts)
  │
  ├─ scripts/run_mafsInRegion_array.sbatch     kent mafsInRegion, per-locus MAF from chrN.maf.gz
  │    (SLURM array 1-109)                     -> maf_output/<part>/<decoyID>.maf
  │
  ├─ scripts/process_maf_spliceai_v2.5.py      msa_view MAF->FASTA (-V on minus strand),
  │    driven by                               SpliceAI 5-model ensemble over each species'
  │    scripts/process_maf_spliceai_array.sbatch  ungapped sequence -> per-part TSV
  │    (SLURM array 1-109)                     -> splice_out/<part>_spliceai_scores.tsv
  │
  └─ scripts/combine_spliceai_tsv.py           concatenate the 109 per-part TSVs
                                               -> results/all_spliceai_scores_v3.1.tsv
```

Then, for the figure:

```
results/all_spliceai_scores_v3.1.tsv + data/tree_to_clade_mapping.tsv
  └─ scripts/hnrnph1decoyposition.Rmd          parse cells, filter decoyID == HNRNPH1_179620582,
                                               plot pos2 (x) vs max2 (y), colour by clade
                                               -> results/hnrnph1_max2_pos2_scatterplot.pdf
                                               -> data/Decoys/hnrnph1_spliceaispecies.tsv
```

## Which scan script produced it — and why it is *not* v3

The filename says `v3.1`, but that is a **dataset/run label, not a script
version**. The scoring script was `process_maf_spliceai_v2.5.py`. Three
independent lines of evidence:

**1. Dates rule v3 out.** `process_maf_spliceai_v3.py` was added in `722cda3` on
**2026-02-25** — 19 days *after* the data file was committed (2026-02-06). The
scripts that existed before the data are `process_maf_spliceai.py` and
`_v2.py` (2026-02-02), `_v2.4.py` and `_v2.5.py` (2026-02-03).

**2. v3 consumes this file, it does not produce it.** Its own docstring:

> *Version 3.0: Uses positional information from v3.1 TSV to extract local 50-nt
> windows and recalculate max1 and max2 scores for each species*

It requires `--v31_tsv` as a mandatory input, and it copies `pos1`, `pos2` and
`d` through unchanged, rewriting only `max1`/`max2`. It is a **re-scoring pass**
layered on top of this table. It cannot be the table's origin.

**3. The `pos1` distribution discriminates v2.4 from v2.5.** v2.4 searches for
`max1` in the first 60 nt (`positions [0, 59]`); v2.5 extends this to 100 nt,
with a length-adjusted branch `int(100 * n/v)` for species longer than hg38. In
the data:

- 457,815 cells carry a valid `pos1`
- **14,387 have `pos1` ≥ 60** — impossible under v2.4's hard cap at 59
- 17 have `pos1` ≥ 100 (max 116) — consistent only with v2.5's `100 * n/v` branch

Reproduce:

```bash
tail -n +2 results/all_spliceai_scores_v3.1.tsv \
  | grep -oE '\[[-0-9.]+\]\[[-0-9]+\]\[' | grep -oE '\]\[[-0-9]+' | tr -d '][' \
  | awk '$1>=60' | wc -l
```

## Status of the v3 re-scoring pass: written, never run in production

`process_maf_spliceai_v3.py` is restored here in full and is byte-exact (see
below), but the evidence says **its output was never used for any committed
result**:

- Its driver `scripts/process_maf_spliceai_array_v3.sbatch` was committed with an
  unfilled placeholder — `V31_TSV="PLACEHOLDER_PATH_TO_V31_TSV_FILE"` above a
  `# TODO: Fill in the path to your v3.1 TSV file`. The script `exit 1`s if that
  path does not resolve, so as committed it could not run.
- It would have written to `splice_out_v3/*_spliceai_scores_v3.tsv`. No such
  files exist anywhere in the repository or its history.
- The leftover `scripts/__pycache__/process_maf_spliceai_v3.cpython-312.pyc`
  embeds `co_filename = intronretention/hnRNPH1_IR_MAF/scripts/process_maf_spliceai_v3.py`
  — a local path, not the HPC `/scratch/.../software/` path the sbatch uses. A
  `.pyc` is written on *import*, not on `python script.py`, so this records local
  interactive/test importing, not a production run.

**Consequence for the manuscript:** the figure's `max1`/`max2` values are v2.5
outputs. Describe the method as the v2.5 whole-intron scan. Do not cite v3 as the
source of these scores.

## Restoration integrity

`process_maf_spliceai_v3.py` was deleted in `aff0a2d` (2026-03-20) and is
restored here from blob `0bf1749ba84f91c8a48861c4a97e28425848e4fd` (commit
`722cda3`). It is verified byte-exact against the surviving `.pyc`:

- `.pyc` header records source size **20659 bytes**; restored file is **20659 bytes**
- `.pyc` header records source mtime **2026-02-24 20:04:13 UTC**
- Recompiling the restored source under CPython 3.12 reproduces **all 19 code
  objects with identical `co_code`, `co_names`, `co_varnames` and constants**

```bash
python3 - <<'EOF'
import marshal
orig = marshal.loads(open("scripts/__pycache__/process_maf_spliceai_v3.cpython-312.pyc",'rb').read()[16:])
new  = compile(open("scripts/process_maf_spliceai_v3.py").read(), orig.co_filename, 'exec')
def sig(c, out):
    out.append((c.co_name, c.co_argcount, c.co_names, c.co_varnames,
                tuple(x for x in c.co_consts if not hasattr(x,'co_code')), c.co_code))
    for k in c.co_consts:
        if hasattr(k,'co_code'): sig(k, out)
    return out
print("EXACT BYTECODE MATCH:", sig(orig,[]) == sig(new,[]))
EOF
```

This is an exact restoration from git, not a reimplementation.

> **Do not import or byte-compile the restored scripts in place.** `import`,
> `python -m py_compile`, and `python -m compileall` all rewrite
> `scripts/__pycache__/*.pyc`, which would destroy the Feb-2026 artefact the check
> above depends on. The `.pyc` is tracked, so an accidental overwrite is
> recoverable with `git checkout -- scripts/__pycache__/`, but check `git status`
> before committing. The snippet above is safe: it reads the `.pyc` and calls
> `compile()` in memory, and writes nothing.

## Files restored 2026-08-07

All recovered from git history at the commits noted; none were rewritten.

| File | From |
|---|---|
| `scripts/split_chromosomes.sh` | `722cda3` |
| `scripts/run_mafsInRegion_array.sbatch` | `da16920` |
| `scripts/run_ctrl_mafsInRegion_array.sbatch` | `da16920` |
| `scripts/process_maf_spliceai.py` | `74762a5` |
| `scripts/process_maf_spliceai_v2.py` | `74762a5` |
| `scripts/process_maf_spliceai_v2.4.py` | `74762a5` |
| `scripts/process_maf_spliceai_v2.5.py` | `74762a5` ← **produced the data** |
| `scripts/process_maf_spliceai_array.sbatch` | `da16920` |
| `scripts/combine_spliceai_tsv.py` | `da16920` |
| `scripts/combine_spliceai_v2.3.py` | `da16920` |
| `scripts/process_maf_spliceai_v3.py` | `722cda3` (re-scoring pass, unused) |
| `scripts/process_maf_spliceai_array_v3.sbatch` | `722cda3` (placeholder, never run) |
| `scripts/stream_maf_region.py` | `722cda3` |
| `scripts/test_stream_maf_hnrnph1.sh` | `722cda3` |

## Known provenance gaps

Stated explicitly so they are not silently papered over:

1. **The committed `process_maf_spliceai_array.sbatch` points at
   `process_maf_spliceai.py`, not `_v2.5.py`.** Its `PYTHON_SCRIPT` variable was
   never updated across the whole of its git history (`da16920` → deletion in
   `e999fee`). The production runs used HPC-local copies under
   `/scratch/prj/ppn_rnp_networks/users/mike.jones/software/`. The v2.5
   attribution above rests on the `pos1` evidence, not on the sbatch.
2. **No SLURM logs are committed**, so exact run dates, wall times and the
   SpliceAI/TensorFlow package versions of the production run are not recoverable
   from this repository.
3. **The upstream MAF and the input BED are not in the repository** (they live on
   HPC scratch). `data/HNRNPH1_179620582.maf` and `_470way.maf` are committed but
   are 0 bytes.
