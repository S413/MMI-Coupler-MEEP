# MMI Coupler MEEP

## Overview

This project focuses on simulating a 2D model of a Multi-Mode Interference (MMI) coupler based on a 3D design from a scientific paper. The simulation is conducted using MEEP (MIT Electromagnetic Equation Propagation), a popular open-source software for simulating electromagnetic systems.

## Table of Contents

- [Overview](#overview)
- [Background](#background)
- [Usage](#usage)
- [Acknowledgments](#acknowledgments)

## Background

MMI couplers are essential components in photonic integrated circuits, enabling the manipulation of light signals in various applications. This project aims to replicate the findings and design presented in a specific scientific paper by conducting 2D simulations. The reduced dimensionality helps in simplifying the computational requirements while still providing insightful results.

## Usage

To run the simulation, follow these steps:

1. Open a terminal and navigate to the project directory.
2. Run the main simulation script and plot + animation:
    ```sh
    python meep_taper_vertices.py tapered_waveguide --plot --gif
    ```

## Acknowledgments

This project is based on the design and research presented in the following scientific paper:
- [Title of the Scientific Paper]
- [Authors]
- [Journal Name, Year, DOI]

---
