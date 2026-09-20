# Contributing to GBS-Genomics-Pipeline

Thank you for your interest in contributing. This project follows open-source best practices for reproducible bioinformatics software.

## How to Contribute

1. **Fork** the repository on GitHub
2. **Create a branch** for your feature or fix: `git checkout -b feature/my-improvement`
3. **Make changes** following existing code conventions
4. **Test** your changes:
   ```bash
   nextflow config
   nextflow run pipeline.nf -stub-run -profile test --curated_dir test_data
   ```
5. **Commit** with a clear message describing the *why*
6. **Push** and open a Pull Request

## Code Standards

- Preserve the existing DSL2 modular architecture
- Do not hard-code paths, usernames, or institution-specific settings
- Use configurable `params` for all user-facing options
- Add informative logging to new processes
- Document new parameters in `docs/PARAMETERS.md`
- Never fabricate biological results in visualization or test outputs

## Reporting Issues

Include:
- Nextflow version (`nextflow -version`)
- Profile used (`docker`, `conda`, etc.)
- Full command line
- Relevant log excerpts from `results/Logs/`
- Input data description (not the data itself)

## Module Development

New modules should:
- Live in `modules/` with a descriptive `.nf` filename
- Include a header comment describing inputs, outputs, and purpose
- Use consistent `publishDir` under the standard output structure
- Declare appropriate `label` for resource allocation

## Scientific Integrity

This pipeline is intended for publication-quality research. Contributors must ensure:
- All annotations come from validated tool outputs
- Missing data produces explicit skip reports, not synthetic placeholders
- Version pins are updated in both Conda envs and Dockerfile when tools change
