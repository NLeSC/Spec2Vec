# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]


## [0.3.0] - 2020-05-13
### Added
- Spectrum, Scores class, save_to_mgf, load_from_mgf, normalize_intensities, calculate_scores https://github.com/matchms/matchms/pull/66 and https://github.com/matchms/matchms/pull/67 https://github.com/matchms/matchms/pull/103 https://github.com/matchms/matchms/pull/108 https://github.com/matchms/matchms/pull/113 https://github.com/matchms/matchms/pull/115 https://github.com/matchms/matchms/pull/151 https://github.com/matchms/matchms/pull/152 https://github.com/matchms/matchms/pull/121 https://github.com/matchms/matchms/pull/154 https://github.com/matchms/matchms/pull/134 https://github.com/matchms/matchms/pull/159 https://github.com/matchms/matchms/pull/161 https://github.com/matchms/matchms/pull/198
- Spikes class https://github.com/matchms/matchms/pull/150 https://github.com/matchms/matchms/pull/167
- Anaconda package https://github.com/matchms/matchms/pull/70 https://github.com/matchms/matchms/pull/68 https://github.com/matchms/matchms/pull/181
- Sonarcloud https://github.com/matchms/matchms/pull/80 https://github.com/matchms/matchms/pull/79 https://github.com/matchms/matchms/pull/149 https://github.com/matchms/matchms/pull/169
- Normalization filter https://github.com/matchms/matchms/pull/83
- SpeciesString filter https://github.com/matchms/matchms/pull/181
- Select by relative intensity filter https://github.com/matchms/matchms/pull/98
- Select-by capability based on mz and intensity https://github.com/matchms/matchms/pull/87
- Default filters https://github.com/matchms/matchms/pull/97
- integration test https://github.com/matchms/matchms/pull/89 https://github.com/matchms/matchms/pull/147 https://github.com/matchms/matchms/pull/156 https://github.com/matchms/matchms/pull/194
- cosine greedy similarity function https://github.com/matchms/matchms/pull/112
- parent mass filter https://github.com/matchms/matchms/pull/116 https://github.com/matchms/matchms/pull/122 https://github.com/matchms/matchms/pull/158
- require_minimum_number_of_peaks filter https://github.com/matchms/matchms/pull/131 https://github.com/matchms/matchms/pull/155
- reduce_to_number_of_peaks filter https://github.com/matchms/matchms/pull/209
- inchi filters https://github.com/matchms/matchms/pull/145 https://github.com/matchms/matchms/pull/127 https://github.com/matchms/matchms/pull/181
- losses https://github.com/matchms/matchms/pull/160
- vesion string checks https://github.com/matchms/matchms/pull/185
- Spec2Vec https://github.com/matchms/matchms/pull/183 https://github.com/matchms/matchms/pull/165 
- functions to verify inchies https://github.com/matchms/matchms/pull/181 https://github.com/matchms/matchms/pull/180
- documentation using radthedocs https://github.com/matchms/matchms/pull/196 https://github.com/matchms/matchms/pull/197
- build status badges https://github.com/matchms/matchms/pull/174
- vectorize spec2vec https://github.com/matchms/matchms/pull/206

### Changed
- Seperate filters https://github.com/matchms/matchms/pull/97
- Translate filter steps to new structure (interpret charge and ionmode) https://github.com/matchms/matchms/pull/73
- filters returning a new spectrum https://github.com/matchms/matchms/pull/100
- Flowchart diagram https://github.com/matchms/matchms/pull/135
- numpy usage https://github.com/matchms/matchms/pull/191
- consistency of the import statements https://github.com/matchms/matchms/pull/189

### Fixed
-

### Removed
-


## [0.2.0] - 2020-04-03
### Added
- Anaconda actions


## [0.1.0] - 2020-03-19
### Added
- This is the initial version of Spec2Vec from https://github.com/iomega/Spec2Vec


[Unreleased]: https://github.com/olivierlacan/keep-a-changelog/compare/v0.3.0...HEAD

[0.3.0]: https://github.com/olivierlacan/keep-a-changelog/compare/v0.2.0...v0.3.0

[0.2.0]: https://github.com/olivierlacan/keep-a-changelog/compare/v0.1.0...v0.2.0

[0.1.0]: https://github.com/olivierlacan/keep-a-changelog/releases/tag/v0.1.0