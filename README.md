# Sizeb_NPZD - Size Compositions of phytoplankton community
A differential equation-based model to study the changes in the size structure of lake phytoplankton communities. The model is embedded with multiple allometric relationships that produce an ecological trade-off favouring small and large phytoplankton at different conditions.



## Model description
The size-based model is succeeded from the well-established Nutrient-Phytoplankton-Zooplankton-Detritus (NPZD) framework (_sensu_ Fasham et al., 1990). This model consist of different size classes of phytoplamkton (Pi) that are subject to grazing from different size classes of zooplankton (Zj). The phytoplankton growth is limited by light and nutrient and is dependent on temperature. 

The model focuses on capture size-dependent bottom-up and top-down interactions through allometric scaling relationships of phytoplankton growth and zooaplankton grazing.



## How to run the model
Running the model requires two scripts in the 'Model run' folder, the 'SizebNPZD_v0.py' and the 'ModRun.py'. The first script is for decribing the model while the second script is for running the model.

Other than these two codes, one would also need the data for environmental or physical forcing in the model, namely the temperature, the irradiance and the nutrient concentration throughout the year.

Default/example forcing data for this model can be found in the 'Data' folder.



## How to generate and analysize results
The results produced from the model is saved as a multi-dimensional arrays into NetCDF format. 



## Related publications
For more detailed descriptions or an example of the model please refer to the related publication of this model (under preparation).
