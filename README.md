# Sizeb_NPZD - Size Compositions of phytoplankton community
A differential equation-based model to study the changes in the size structure of lake phytoplankton communities. The model is embedded with multiple allometric relationships that produce an ecological trade-off favouring small and large phytoplankton at different conditions.

The allometric relationships considered for four eco-physiological rates in the model are:

$$\mu_{max}(S_i^P) = \beta_{\mu_{max}}\cdot (S_i^P)^{\alpha_{\mu_{max}}}$$

$$K_n(S_i^P) = \beta_{K_n}\cdot (S_i^P)^{\alpha_{K_n}}$$

$$I_{max}(S_j^Z) = \beta_{I_{max}}\cdot (S_j^Z)^{\alpha_{I_{max}}}$$

$$P_{opt}(S_i^P, S_j^Z) = \beta_{P_{opt}}\cdot (S_j^Z)^{\alpha_{P_{opt}}}$$

representing maximum growth rate, $\mu_{max}(S_i^P)$, and half-saturation for nutrient uptake, $K_n(S_i^P)$, for phytoplankton size class $i$, and maximum ingestion rate, $I_{max}(S_j^Z)$, and optimal prey size, $P_{opt}(S_i^P, S_j^Z)$, for zooplankton size class $j$ respectively.

## Model description
The size-based model is adapted from the well-established Nutrient-Phytoplankton-Zooplankton-Detritus (NPZD) framework (_sensu_ Fasham et al., 1990 and Armstrong, 1994). This model consists of different size classes of phytoplankton ($P_i$) who are subject to grazing by two size classes of zooplankton ($Z_1$, $Z_2$). The phytoplankton growth is limited by light and nutrient and is dependent on temperature. 

The model focuses on capturing size-dependent bottom-up and top-down interactions through allometric scaling relationships of phytoplankton growth and zooaplankton grazing.

<p align="center">
  <img width="560" alt="Figure1_v3" src="https://github.com/Debbcwing/Sizeb_NPZD/assets/51200142/a853f676-68e2-45b4-99e0-8c42f0c4c1a6">

## How to run the model
Running the model requires two scripts in the 'Model run' folder, the 'SizebNPZD_v0.py' and the 'ModRun.py'. The first script is for decribing the model while the second script is for running the model.

Other than these two scripts, one would also need datasets for environmental/physical forcing to the model, namely, the temperature, the irradiance and the nutrient concentration throughout the year. Default/example forcing data for this model can be found in the 'Data' folder.



## How to generate and analysize results
The results produced from the model are saved as a multi-dimensional array in NetCDF format and archived in Zenodo (doi: 10.5281/zenodo.7431914)


## Related publications
For more detailed model descriptions, formulations, or an application example of the model, please refer to the related publication of this model:

To, S., E. Acevedo‐Trejos, S. Chakraborty, F. Pomati, and A. Merico. 2024. Grazing strategies determine the size composition of phytoplankton in eutrophic lakes. _Limnology & Oceanography_ 9999:1-10. doi:10.1002/lno.12538

