# Making canonical scripts portable

Several analysis scripts were originally written against the local
`Benchmark_Program` directory. For the public release, **do not change scientific logic** merely to make
paths portable.

Preferred approaches, in order:

1. Add a command-line `--root` argument if the script already supports one.
2. Derive project root relative to the script location when this is already part of the canonical design.
3. Read paths from `config/paths.example.yaml` or environment variables.
4. If a script must remain exactly frozen for provenance, keep the frozen file unchanged and add a small
   wrapper script that supplies portable paths.

Always preserve the frozen canonical file and its SHA256. If a portability patch is made, version it as a
separate release script and document the diff.
