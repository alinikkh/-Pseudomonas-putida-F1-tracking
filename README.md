# Pseudomonas-putida-F1-tracking
Pseudomonas-putida-F1-Tracking

This repository contains the tracking code for analyzing the motility of Pseudomonas putida F1 cells. The code is designed to process experimental trajectories and detect bacterial runs and tumbles, facilitating quantitative analysis of bacterial motion under various experimental conditions.

Input Requirements

The code requires an Excel file with three columns:

Frame number – the sequential frame index from your recorded video.

X-position – horizontal position of the cell in pixels.

Y-position – vertical position of the cell in pixels.

To convert the data into physical units, users must provide:

Microscope resolution: to convert pixel coordinates into micrometers.

Frame rate: to convert frame indices into seconds.

Tumble Detection

The algorithm identifies bacterial tumbling events based on criteria described in our paper. These criteria can be adjusted and visually verified to ensure accurate detection for your specific dataset.

Usage

Prepare your trajectory data in the required Excel format.

Input your microscope resolution and frame rate in the code.

Run the script to analyze trajectories and extract run and tumble events.

Visualize and validate the results to ensure accuracy.

Notes

The code has been tested on trajectories obtained from videos recorded at 40 fps.

Users can modify the tumble detection parameters if the experimental setup or bacterial strain differs.

Output includes processed trajectories and a summary of run/tumble statistics.

Citation
    
If you use this code for your research, please cite:
Doan et al., “Stabilizing by steering: Enhancing bacterial motility by non-uniform diffusiophoresis,” Supporting Information.
