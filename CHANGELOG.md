# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),  
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [3.1.2] - 2025-06-03

### Added
- CLI (`tobacco`) with multiple subcommands (`--make_MOF`, `--gen_geometries`, etc.).
- Functionality to download topology database using `--get_topols_db`.
- New geometry generator from atomic point group and metal.
- Support for `.com` (Gaussian) to SBU conversion.
- `--build_sbu node/edge`, `--make_XX_edge`, and other SBU tools.
- Compatibility with Python 3.8+ using virtual environments.

### Changed
- Refactored main entry point from `tobacco.py` to CLI tool.
- Internal reorganization of files and modules.
- Improved modularity for MOF and geometry building scripts.

### Deprecated
- Direct execution of `tobacco.py` is discouraged in favor of the CLI interface.

### Removed
- Legacy files and examples not compatible with new structure.
- Old installation notes based on Python 2.7 (now deprecated).

### Fixed
- Bug in topology parser for some edge cases.
- CIF output atom order consistency.

---

## [3.0.0] - 2019

Initial public release by Ryther Anderson, Yamil Colón, and Diego Gómez-Gualdrón.