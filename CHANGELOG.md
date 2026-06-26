# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html)
with an additional version indicator for modifications of this fork.

## [Unreleased]

### Added <!-- New features -->

-

### Changed <!-- Changes in existing functionality -->

-

### Deprecated <!-- Soon-to-be removed features -->

-

### Removed <!-- Now removed features -->

-

### Fixed <!-- Any bug fixes -->

-

### Security <!-- In case of vulnerabilities -->

-

## [v3.5.1.0]

Release Date: March, 2025

This release of the LaSRC code has been modified by the NASA-IMPACT team as part of
maintenance for the HLS project based off of the USGS's version `v3.5.1`.
We are indicating this change by appending a "maintainer version" to the version tag
(`[major].[minor].[patch].[maintainer]`). This is done to prevent potential confusion
with releases from the upstream.

### Fixed

- Support Sentinel-2C and -2D by updating hard coded platform values
- Include `PROC_ALL_BANDS` definition
- Mask a pixel if _any_ band has invalid data [#17](https://github.com/NASA-IMPACT/espa-surface-reflectance/pull/17)


[Unreleased]: https://github.com/NASA-IMPACT/espa-surface-reflectance/tree/v3.5.1.0...HEAD
[v3.5.1.0]: https://github.com/NASA-IMPACT/espa-surface-reflectance/tree/v3.5.1.0
