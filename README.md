# Robust Profile Likelihood (RPG) for Robust Optimal Design

This repository contains the code and data used to generate **Robustness Profile Graphics (RPGs)** for analyzing optimal and robust experimental designs. The RPG provides insight into the contribution of individual design points to overall design efficiency, allowing for a more detailed assessment of design robustness under potential data loss.

The study extends previous methods by using **particle swarm optimization** to generate robust designs, and visualizes them at both the overall and observation levels.

---

## Repository Structure

Study_Collection/
├── Study_Code.R # Main R script to generate plots and RPGs
├── K=1/ # Design scenarios for K=1
│ ├── N=4/ ... N=9/ # Subfolders containing CSVs for 500 iterations (G and I robust criteria)
│ ├── best_designs.csv # Contains best designs for I, G, robust I, robust G, and Central Composite Design
│ └── RPGs.pdf # PDFs showing all six RPGs
├── K=2/ # Same structure for K=2
├── K=3/ # Same structure for K=3
