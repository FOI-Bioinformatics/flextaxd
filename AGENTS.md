# Project notes

- Read `FLEXTAXD.md` (project root) before reviewing or fixing anything — it describes the intended functionality of FlexTaxD and what must be preserved (NCBI taxids as gold standard, name-based merging, GTDB swap, clear_database intent).
- Some code looks wrong but is intentional: e.g. `split("\t|\t")` in ReadTaxonomyNCBI is the NCBI dmp field separator, and the `database=` parameters in DatabaseConnection methods are unused (passing them is harmless). Check FLEXTAXD.md and surrounding comments before treating something as a bug.
- Issues live in `ISSUES.md`, one line per issue, managed by `issues.sh` (don't reorder or renumber by hand).
- The git repository is this folder (`flextaxd/`), the folder above it is not a repo. Write the suggested commit message to `flextaxd/.commit-msg` here, not to the folder above — the user commits it with `opencode-sandbox commit` from here.
