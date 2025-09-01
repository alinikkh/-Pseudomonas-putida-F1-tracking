# Pseudomonas-putida-F1-tracking
# 🦠 Pseudomonas-putida-F1-Tracking

[![MATLAB](https://img.shields.io/badge/MATLAB-R2020b-blue)](https://www.mathworks.com/products/matlab.html)
[![License](https://img.shields.io/badge/License-MIT-green)](LICENSE)
[![GitHub Repo Size](https://img.shields.io/github/repo-size/yourusername/Pseudomonas-putida-F1-Tracking)](https://github.com/yourusername/Pseudomonas-putida-F1-Tracking)

This repository contains the **MATLAB code** for analyzing the motility of *Pseudomonas putida* F1 cells. The code detects bacterial **runs and tumbles** from experimental trajectories and enables quantitative analysis of bacterial motion under various experimental conditions.
## 📚 Update
  ▶ Python code has been developed and is under test. (Will be added soon).
  ▶ The New code will give more organized figures.
  ▶ All Fig files will be added to a separate folder.


---

## ⚡ Features

- ✅ Analyze bacterial trajectories from microscopy videos  
- ✅ Detect runs and tumbles with adjustable criteria  
- ✅ Convert pixel coordinates to physical units using microscope resolution  
- ✅ Visualize trajectories and generate color-coded straightness maps  
- ✅ Generate movies showing bacterial motion under different conditions  

---

## 📂 Input Requirements

The code requires an **Excel file** with **three columns**:  

| Column | Description |
|--------|-------------|
| 1      | Frame number (sequential frame index) |
| 2      | X-position (in pixels) |
| 3      | Y-position (in pixels) |

**Conversion requirements:**  
- **Microscope resolution** → convert pixels to micrometers  
- **Frame rate** → convert frame numbers to seconds  

---

## 🔍 Tumble Detection

- Based on criteria described in our paper  
- Parameters can be **adjusted** and **visually verified** for each dataset  
- Ensures accurate detection of bacterial runs and tumbles  

---

## 🚀 Usage

1. Open MATLAB and navigate to the repository folder  
2. Load your trajectory Excel file  
3. Set microscope resolution and frame rate in the code  
4. Run the main script to analyze trajectories  
5. Visualize and validate results  

---

## 🎥 Visualization

The code can generate clear visual outputs:  

- **Run trajectories** – displacement and orientation of each bacterium  
- **Color-coded straightness** – highlights directional motion  
- **Movies** – animations of bacterial motion under different conditions  

<img width="3020" height="1457" alt="runtimedistribution" src="https://github.com/user-attachments/assets/d2aaa30b-762b-44c4-8196-27cfbb00e7ec" />


**Example outputs:**     
- **Movie S4–S5:** Run trajectories under uniform and gradient conditions, experimental and simulated  
---

## 📝 Notes

- Tested with **MATLAB R2024b**   
- Tumble detection parameters can be modified for different experimental setups and bacteria  
- Output includes processed trajectories, run/tumble events, and visualizations  

---

## 📌 Citation

If you use this code in your research, please cite our paper:  
*Doan et al., “Stabilizing by steering: Enhancing bacterial motility by non-uniform diffusiophoresis,” Supporting Information.*

---

## ⚙️ Requirements

- MATLAB R2024b or later  
- MATLAB toolboxes:  
  - `Statistics and Machine Learning Toolbox` (optional for advanced analysis)  
  - `Image Processing Toolbox like Image J` (if working with videos)  

