# Co-ordinated Beamforming for Interference Management in CoMP-Based Integrated Terrestrial and Non-Terrestrial Network

This repository contains MATLAB simulations for a final year project at Kyambogo University by Arnold Odongo and Hafswa Ndagire Babirye, submitted in June 2025 under the supervision of Dr. Denise J. Birabwa. The project analyzes the Signal-to-Interference-plus-Noise Ratio (SINR) performance of Coordinated Multi-Point (CoMP) and Non-CoMP beamforming techniques in an Integrated Terrestrial and Non-Terrestrial Network (ITNTN). The simulations compare Regularized Zero-Forcing (RZF), Zero-Forcing (ZF), and Maximum Ratio Transmission (MRT) precoding methods under varying user counts..

## Project Overview

The project models a two-tier Random Access Network (RAN) with 4 base stations (3 Macro Base Stations (MBS) and 1 Low-Altitude Platform Base Station (LAP-BS)) in a 5 km x 5 km area. Users move according to a random walk mobility model. The simulations evaluate Complementary Cumulative Distributon Function(CCDF) of SINR as well as median SINR as a function of the number of users (10, 20, 30, 40, 50) over 1000 time slots, considering path loss, shadowing, and fading (Rayleigh for MBS, Rician for LAP-BS). The work addresses interference management challenges in ITNTNs, focusing on co-tier and cross-tier interference.

## Files

- **`ComparingCoMPRZFvsNonCoMPRZF.m`**:
  - Compares RZF beamforming for both CoMP and non-CoMP scenarios.
  - Simulates SINR for 50 users over 1000 time slots.
  - Outputs a CCDF plot of SINR (dB) for CoMP and Non-CoMP RZF.


- **`MedianSINRComparisonCoMPvsNon-CoMP.m`**:
  - Compares CoMP RZF and Non-CoMP RZF beamforming for varying user counts (10, 20, 30, 40, 50).
  - Outputs a plot of median SINR (dB) vs. number of users.


- **`ComparingRZFZFandMRTPrecodingwithinCoMP.m`**:
  - Compares RZF, ZF, and MRT precoding within CoMP for a fixed 50 users.
  - Generates a network topology plot and SINR CCDF plot.

  
- **`MedianSINRComparisonRZFZFandMRT.m`**:
  - Compares RZF, ZF, and MRT precoding within CoMP for varying user counts (10 to 50).
  - Outputs a plot of median SINR (dB) vs. number of users.
  

## Prerequisites

- **MATLAB**: Version R2024a or later recommended.
- **System Requirements**: Standard desktop/laptop with sufficient memory for simulations (1000 time slots per user count).
