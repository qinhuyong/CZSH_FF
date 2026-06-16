# Package Manifest

| File/folder | Purpose | Executable | Generated output | Notes |
| --- | --- | --- | --- | --- |
| `README.md` | Package overview and scope | no | no | Start here for package purpose. |
| `VERSION.txt` | Version and source tag metadata | no | no | Records v1.11 package identity. |
| `MANIFEST.md` | Package content map | no | no | This file. |
| `CITATION_NOTE.md` | Citation placeholders and citation reminders | no | no | Does not invent reference formats. |
| `LICENSE_NOTE.md` | License status note | no | no | License should be finalized before archival release. |
| `QUICK_START.md` | Minimal command sequence | no | no | Uses package-relative paths. |
| `REPRODUCIBILITY_COMMANDS.md` | Lightweight verification commands | no | no | Includes expected `07` target list. |
| `pyCSH_Zn/` | Workflow source snapshot and documentation | mixed | no | Excludes `output_Y/`, caches, and logs. |
| `pyCSH_Zn/examples/` | Workflow entry-point scripts | yes | no | Python scripts are unchanged snapshot files. |
| `pyCSH_Zn/forcefields/` | CementFF4-Zn parameter database and builders | mixed | no | Force-field parameters are unchanged. |
| `pyCSH_Zn/lammps_templates/` | LAMMPS input template generation | mixed | no | Used by workflow scripts. |
| `pyCSH_Zn/postprocess/` | RDF, coordination, angle, and contact analysis | mixed | no | Optional post-processing helpers. |
| `pyCSH_Zn/FINAL_WORKFLOW.md` | Full workflow reference | no | no | Detailed command and label tables. |
| `pyCSH_Zn/MANUSCRIPT_METHODS_DRAFT.md` | Methods draft | no | no | Manuscript text starting point. |
| `pyCSH_Zn/MANUSCRIPT_FIGURE_PLAN.md` | Figure plan | no | no | Figure concepts and data sources. |
| `pyCSH_Zn/SUPPLEMENTARY_WORKFLOW_TABLES.md` | Supplement tables | no | no | Script, label, and output column tables. |
| `pyCSH_Zn/RELEASE_AUDIT_v1.9.md` | v1.9 audit | no | no | Documents frozen boundaries. |
| `checks/PAPER_CODE_AUDIT.md` | Package audit | no | no | Confirms exclusions and boundaries. |
