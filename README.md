## Landscape of Epithelial Mesenchymal Plasticity as an emergent property of coordinated teams in regulatory networks

This resource provides the R code required for all the data simulations and analysis done in [**Landscape of Epithelial Mesenchymal Plasticity as an emergent property of coordinated teams in regulatory networks**](https://www.biorxiv.org/content/10.1101/2021.12.12.472090v2). The raw data and processed data generated in the analysis can be found at <https://drive.google.com/drive/folders/1bkVEQ7Wn4shB-kd8sfW_WKs-XeNZs_Vu?usp=sharing>. The folder structure can be understood by studying the setupScript.R file. 

The simulation tools used are as follows:

* **Discrete Modelling**: Asynchonous Update on Ising model of WT, perturbed and random EMP networks. This is done by the set of codes given [here](https://github.com/askhari139/Boolean.jl). Keeping with the traditional Ising model, states are represented as -1 and 1.
* **Continuous Modelling**: **RACIPE** is used to generate an ensemble of continuous models of WT-EMP networks. For this we have used a modified version of [**RACIPE-1.0**](https://github.com/simonhb1990/RACIPE-1.0) package which can be found [here](https://github.com/csbBSSE/Gene_Network_Modelling/releases/download/v2.28/Multithreaded_Racipe_2.28.zip).

A brief description of the analysis:

**1.** Generate random networks (single edge deletions and multi-edge deletions) for the WT EMP topo files

**2.** Simulate all the WT networks (and perturbed networks for Boolean) using Boolean and Continuous formalisms listed above. The boolean simulations result in the steady state frequencies and frustration. For RACIPE data, frustration would have to be calculated.

**3.** Calculate the group strengths of all the networks.

**4.** Generate correlation matrices for all the network simulations, both RACIPE and Boolean.

**5.** Calculate coherence of all steady states, as well as multi-node perturbation coherence.

**6.** Calculate State strength, EMT score for all states of each network and attach phenotype labels for each of the states. 

*** The files generated using these codes are available on the google drive [here](https://drive.google.com/drive/folders/1bkVEQ7Wn4shB-kd8sfW_WKs-XeNZs_Vu?usp=sharing).

<!---
### Figures
All the figures presented in the paper including Supplementary Figures (apart from the Schematics and Network Representation) are provided in the [``Figures``](https://github.com/csbBSSE/CSB-SCLC/tree/master/Figures) folder. Details of reproducing the figures are briefly given in the in the ``README`` file in each folder containing the subfigures.

### Simulation Data
**Boolean Simulation data** generated using [``Fast-Bool``](https://github.com/csbBSSE/CSB-SCLC/tree/master/Additional_Codes/Fast-Bool) and **Edge Perturbation data** generated using [``Edge_Perturbation``](https://github.com/csbBSSE/CSB-SCLC/tree/master/Additional_Codes/Edge_Perturbation) are provided in the **Simulation_Data** folder. Since **RACIPE** simulation data files are quite huge, they are uploaded to this [drive link](https://drive.google.com/drive/folders/1PKs5vHkXCoJm9Wcg7P4nBPdPrFJCxJ5B?usp=sharing).

### Additonal Codes
This folder contains all the codes required for data analysis and simulation of WT-SCLC network. These are the basic framework of codes which some scripts used for figure production relies on.

### How to reproduce the plots?
**1.** Clone the GitHub Repository
```
git clone https://github.com/csbBSSE/CSB-SCLC
```
**2.** Set the working directory to ``CSB-SCLC``

**3.** Install all the Required Python Packages (**Conda** is preferable)
```
while read requirement; do conda install --yes $requirement || pip install $requirement; done < requirements.txt
```
**4.** Go to the folder corresponding to the figure that needs to be reproduced and follow the instructions given there

**5.** Voila! you are done

### Notes
* Some of the codes have an option of running processes in Parallel. Just make sure that you don't give spawn more processes than your CPUs can handle.
* Installing Python packages using **Conda** would be preferable. Using **Intel Python Distribution** gives significant speed boosts in some of the codes.
* Codes like **UMAP_analysis** and **Bool.py** may take longer times. Just be patient and don't **Ctrl+C** it even if you have to wait for some time (_Just Don't do it. Time is precious_)
-->

### Requirements
R(Tested on Version 4.1.2)
Julia(Tested on Version 1.6.2 & 1.4)




 
 
