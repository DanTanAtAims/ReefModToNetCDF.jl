# Converting Reefmod Matfiles to a Singular NetCDF file

## Usage

The file provides two top level functions for conversion.

```julia
# All extracted variables will be written to a single NetCDF file
to_NetCDF_combined(
    "<path/to/matfile/directory>", "<climate scenario>"; out_path="<file/write/path>"
)
```

```julia
# Input and output variables relatife to ADRIA will be written to seperate files
to_NetCDF_combined(
    "<path/to/matfile/directory>", "<climate scenario>"; 
    input_write_path="<file/write/path>", output_write_path="<file/write/path>"
)
```

## Expected Directory Structure

The conversion function assumes all matfiles relating to the given climate scenario are in
the same directory and have the climate scenario contained in there filename.

For example 

```bash
├───reefmod_matfiles
│       sR0_GBR.7.0_herit0.3_SSP119_CNRM-ESM2-1.mat
│       sR0_GBR.7.0_herit0.3_SSP119_EC-Earth3-Veg.mat
│       sR0_GBR.7.0_herit0.3_SSP119_IPSL-CM6A-LR.mat
│       sR0_GBR.7.0_herit0.3_SSP119_MIROC-ES2L.mat
│       sR0_GBR.7.0_herit0.3_SSP119_MIROC6.mat
│       sR0_GBR.7.0_herit0.3_SSP119_MRI-ESM2-0.mat
│       sR0_GBR.7.0_herit0.3_SSP119_UKESM1-0-LL.mat
```

Noting the "SSP119" is the climate scenario and the succeeding text is the earth system
model used.

## Expected Matfile structure

The following keys are assumed to be within the matfile

```julia
"record_applied_DHWs",
"record_applied_cyclones",
"record_applied_bleaching_mortality",
"coral_HT_mean",
"coral_HT_var",
"COTS_juv_densities",
"COTS_adult_densities",
"coral_larval_supply",
"coral_cover_per_taxa",
"nb_coral_offspring",
"nb_coral_adol",
"COTS_larval_supply",
"nb_coral_adult",
"nb_coral_juv",
"nb_coral_recruit",
"reef_shelter_volume_relative"
```

## YAXArray Properties

Useful metadata is also collected from the matfiles and stored in the NetCDF to provide
context and information to users. They can be accessed as follows

```julia
julia> using NetCDF, YAXArrays
# Load dataset
julia> ds = open_dataset("path/to/netcdf.nc")

julia> ds.properties
Dict{String, Any} with 10 entries:
  "climate scenario (RCP)"        => "45"
  "given filename"                => ["sR0_GBR.7.0_herit0.3_SSP245_CNRM-ESM2-1",...] 
  "scenarios per earth model"     => [20, 20, 20, 20, 20, 20, 20, 20, 20, 20]
  "scenario start index"          => [1, 21, 41, 61, 81, 101, 121, 141, 161, 181]
  "scenario end index"            => [21, 41, 61, 81, 101, 121, 141, 161, 181, 201]
  "earth system model order"      => ["CNRM-ESM2-1",]
  "metadata filename"             => ["R0_GBR.7.0_herit0.3_SSP245_CNRM-ESM2-1",...]
  "reef model"                    => "ReefMod Matlab"
  "given filename description"    => "filename of matfile used to create this netcdf"
  "metadata filename description" => "filename contained in the metadata of original matfile"
```
