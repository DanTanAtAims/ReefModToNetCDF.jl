using FilePaths
using NetCDF
using YAXArrays
using MAT
using ProgressMeter

"""
    _get_matfiles(directory::String, scenario::String)::Vector{String}

Get all matfiles relating to the given climate scenario.
"""
function _get_matfiles(directory::String, scenario::String)::Vector{String}
    pos_files = filter(isfile, readdir(directory, join=true))
    pos_files = filter(x -> occursin(".mat", x), pos_files)
    return filter(x -> occursin(scenario, x), pos_files)
end

"""
    _extract_n_scenarios(dataset::Vector{Dict{String, Any}})::Vector{Int}

Extract the number of scenarios run per earth system model. Typically, all will have the same size.
Assume the dimensions of DHWs are [scenario ⋅ locations ⋅ timestep ⋅ [taxonamy]]
"""
function _extract_n_scenarios(dataset::Vector{Dict{String, Any}})::Vector{Int}
    # ADRIA will always require applied DHWs so we use the shape of DHW matrix
    return [size(ds["record_applied_DHWs"])[1] for ds in dataset]
end

"""
    extract_filename(filepath::String)::String

Given the path to a file, return the name of the file without the prefix and suffix.
"""
function _extract_filename(filepath::String)::String
    filename = basename(filepath)
    filename, _ = splitext(filename)
    return filename
end

"""
    _extract_earth_model(original_filename::String)::String 

Extract the name of the earth system model from the filename.

Assumes the following format "*SSP_*<Name of model>.mat"
"""
function _extract_earth_model(original_filename::String)::String
    
    found = match(r"SSP\d*_(.*)", original_filename)
    if isnothing(found)
        @info "Earth system model name could not be found in given filename, returning \"unknown\"."
        return "unknown"
    end

    return String(found.captures[1])
end

"""
    _extract_earth_models(dataset::Vector{Dict{String, Any}}, filepaths::Vector{String})::Vector{String}

Given the a list of datasets extract the earth model from the filename contained in the metadata.
If the original filename is missing from metadata fall back to the current filename.
"""
function _extract_earth_models(dataset::Vector{Dict{String, Any}}, filepaths::Vector{String})::Vector{String}
    @assert length(dataset) == length(filepaths) "Given list of data sets and list of filepaths have different lengths of $(length(dataset)) and $(length(filepaths)) respectively"
    original_filenames = [_extract_original_filename(ds) for ds in dataset]
    
    # Using filename/path if original name not found in metadata
    earth_models = [
        _extract_earth_model(extr == "unknown" ? fn : extr ) 
        for (extr, fn) in zip(original_filenames, filepaths)
    ]
    return earth_models
end

"""
    _extract_climate_scenario(original_filename::String)::String

Extract climate scenario file OutputFileName contained in the metadata of the original matfile
"""
function _extract_climate_scenario(original_filename::String)::String

    found = match(r"SSP\d*", original_filename)
    if isnothing(found)
        @info "Climate scenario could not be found in given filename, returning \"unknown\"."
        return "unknown"
    end
    return found.match
end

"""
    _extract_original_filename(dataset::Dict{String, Any})::String

Extract the original filename contained in the matfiles metadata.
"""
function _extract_original_filename(dataset::Dict{String,Any})::String
    if !haskey(dataset, "OPTIONS")
        @info "Given matfile does not have `OPTIONS` key. Listing original file name as \"unknown\"."
        return "unknown"
    elseif !haskey(dataset["OPTIONS"], "OutputFileName")
        @info "Given matfile \"OPTIONS\" section does not contain `OutputFileName` key. Listing original file name as \"unknown\"."
        return "unknown"
    end

    return dataset["OPTIONS"]["OutputFileName"]
end

"""
    _construct_properties(
        dataset::Dict{String, Any}, mat_filename::String, RCP::String
    )::Dict{String, Any}

Construct the properties for the YAX dataset.
"""
function _construct_properties(
    dataset::Vector{Dict{String,Any}}, filenames::Vector{String}, RCP::String
)::Dict{String,Any}

    earth_models = _extract_earth_models(dataset, filenames)
    scenario_per_model = _extract_n_scenarios(dataset)
    scenario_range_beginning = Vector{Int}(undef, length(filenames))
    scenario_range_ending = Vector{Int}(undef, length(filenames))
    reef_ids = _extract_reef_IDs(dataset[1]) # assume duplicate

    n_scenarios::Int = sum(scenario_per_model)
    n_locations::Int = length(reef_ids) 

    begin_index = 1
    for i ∈ 1:length(filenames)
        scenario_range_beginning[i] = begin_index
        scenario_range_ending[i] = begin_index + scenario_per_model[i] - 1
        begin_index += scenario_per_model[i]
    end

    metadata_filenames = [_extract_original_filename(ds) for ds in dataset]
    given_filenames = [_extract_filename(fn) for fn in filenames]

    reefmodel = "ReefMod Matlab"

    return Dict(
        "climate scenario (RCP)" => RCP,
        "reef model" => reefmodel,
        "metadata filename" => metadata_filenames,
        "metadata filename description" => "filename contained in the metadata of original matfile",
        "given filename" => given_filenames,
        "given filename description" => "filename of matfile used to create this netcdf",
        "earth system model order" => earth_models,
        "scenarios per earth model" => scenario_per_model,
        "scenario start index" => scenario_range_beginning,
        "scenario end index" => scenario_range_ending,
        "reef_ids" => reef_ids,
        "n_scenarios" => n_scenarios,
        "n_locations" => n_locations
    )
end


"""
_extract_reef_IDs(mat_dict::Dict{String, Any})::Vetcor{String}

Extract Reef IDs from dataset metadata.
"""
function _extract_reef_IDs(mat_dict::Dict{String,Any})::Vector{String}
    n_locations = 0

    if !haskey(mat_dict, "META")
        @info "Given matfile does not have 'META' key. Returning empty YAXArray"
        return []
    elseif !haskey(mat_dict["META"], "reef_ID")
        @info "Given matfile META dict does not have 'reef Id' key. Returning index as id for reef area"
        n_locations = mat_dict["META"]["nb_reefs"]
        return string.(1:n_locations)
    end
    (n_locations, _) = size(mat_dict["META"]["reef_ID"])

    return string.(dropdims(mat_dict["META"]["reef_ID"], dims=2))
end

"""
    _extract_lon(mat_dict::Dict{String, Any})::YAXArray

Extract the longitude position of target locations in the dataset.
"""
function _extract_lon(mat_dict::Dict{String,Any})::YAXArray
    if !haskey(mat_dict, "META")
        @info "Given matfile does not have 'META' key. Returning empty YAXArray"
        return YAXArray((), [])
    elseif !haskey(mat_dict["META"], "reef_lon")
        @info "Given matfile META dict does not have 'area_habitat' key. Returning empty YAXArray for reef area"
        return YAXArray((), [])
    end
    (n_locations, _) = size(mat_dict["META"]["area_habitat"])

    axlist::Tuple = (
        Dim{:location}(range(1, n_locations)),
    )

    props::Dict{String,Any} = Dict(
        "location dimension" => 1,
        "variable" => "reef latitude",
    )
    return YAXArray(axlist, dropdims(mat_dict["META"]["reef_lon"], dims=2), props)
end

"""
    _extract_lat(mat_dict::Dict{String, Any})::YAXArray

Extract the latitude position of target locations in the dataset.
"""
function _extract_lat(mat_dict::Dict{String,Any})::YAXArray
    if !haskey(mat_dict, "META")
        @info "Given matfile does not have 'META' key. Returning empty YAXArray"
        return YAXArray((), [])
    elseif !haskey(mat_dict["META"], "reef_lat")
        @info "Given matfile META dict does not have 'area_habitat' key. Returning empty YAXArray for reef area"
        return YAXArray((), [])
    end
    (n_locations, _) = size(mat_dict["META"]["area_habitat"])

    axlist::Tuple = (
        Dim{:location}(range(1, n_locations)),
    )

    props_lat::Dict{String,Any} = Dict(
        "location dimension" => 1,
        "variable" => "reef longitude",
    )
    return YAXArray(axlist, dropdims(mat_dict["META"]["reef_lat"], dims=2), props_lat)
end

"""
    _extract_reef_area(mat_dict::Dict{String, Any})::YAXArray

Extract reef area from from dataset metadata.
"""
function _extract_reef_area(mat_dict::Dict{String,Any})::YAXArray
    if !haskey(mat_dict, "META")
        @info "Given matfile does not have 'META' key. Returning empty YAXArray"
        return YAXArray((), [])
    elseif !haskey(mat_dict["META"], "area_habitat")
        @info "Given matfile META dict does not have 'area_habitat' key. Returning empty YAXArray for reef area"
        return YAXArray((), [])
    end

    (n_locations, _) = size(mat_dict["META"]["area_habitat"])

    axlist::Tuple = (
        Dim{:location}(range(1, n_locations)),
    )

    props::Dict{String,Any} = Dict(
        "location dimension" => 1,
        "variable" => "reef area",
        "units" => "km2"
    )

    return YAXArray(axlist, dropdims(mat_dict["META"]["area_habitat"], dims=2), props)
end

"""
    _extract_nongrazable_area(mat_dict::Dict{String, Any})::YAXArray

Extract nongrazable area from the dataset.
"""
function _extract_nongrazable_area(datasets::Vector{Dict{String,Any}})::YAXArray
    if !all(haskey.(datasets, "nongrazable"))
        @info "Given matfile does not have 'nongrazable' key. Return empty YAXArray for nongrazable"
        return YAXArray((), [])
    end
    final_arr = _join_variable_arr([ds["nongrazable"] for ds in datasets])
    (n_scenarios, n_locations) = size(final_arr)

    axlist::Tuple = (
        Dim{:location}(range(1, n_locations)),
        Dim{:scenario}(range(1, n_scenarios)),
    )

    props::Dict{String,Any} = Dict(
        "variable" => "nongrazable area",
        "units" => "percentage of reef area",
    )

    return YAXArray(axlist, permutedims(final_arr, [2, 1]), props)
end

function equals_shape_excluding(dimension::Int, shape_1::Tuple, shape_2::Tuple)::Bool
    if length(shape_1) != length(shape_2)
        return false
    end
    n_dims = length(shape_1)
    equal = true
    for i in 1:n_dims
        if (i == dimension)
            continue
        end
        equal = equal && (shape_1[i] == shape_2[i])
    end
    return equal
end

"""
    _join_variable_arr(variable_arrs::Vector{Array{<:Real, 3}})::Array{<:Real, 3}
    _join_variable_arr(variable_arrs::Vector{Array{<:Real, 4}})::Array{<:Real, 4}

Concatenate variable arrays together along the scenario dimension. Assumes the following 
dimensions.

Dimensions ordering before reordering with standard AIMS formatting.
[scenario ⋅ locations ⋅ timestep]
[scenario ⋅ locations ⋅ timestep ⋅ group]
"""
function _join_variable_arr(variable_arrs::Vector{<:Array{<:Real}})::Array{<:Real}
    shape = size(variable_arrs[1])
    shps = [equals_shape_excluding(1, size(arr), shape) for arr in variable_arrs]
    if !all(shps)
        ArgumentError("Variable arrays given are not the same shape. Possible dataset mismatch")
    end
    
    fin = reduce((x, y) -> cat(x, y, dims=1), variable_arrs)
    return fin
end

"""
    _extract_variable(dataset::Dict{String, Any}, variable::String)::YAXArray

Extract variable from dataset and return YAXArray.
"""
function _extract_variable(
    dataset::Vector{Dict{String,Any}}, variable_name::String, first_year::Float64
)::YAXArray
    if !all(haskey.(dataset, variable_name))
        @info "$(variable) not found in dataset. Returning empty YAXArray."
        return YAXArray((), [])
    end
    variable_arr = [ds[variable_name] for ds in dataset]
    n_dimensions = length(size(variable_arr[1]))

    if n_dimensions == 3 || n_dimensions == 4
        return _extract_variable(variable_arr, variable_name, Int(first_year))
    else
        @warn "Variable with unexpected dimensions given. Returning empty."
        return YAXArray((), [])
    end
end
function _extract_variable(
    variable_arr::Vector{<:Array{<:Real, 3}}, variable_name::String, first_year::Int
)::YAXArray
    final_arr = _join_variable_arr(variable_arr)
    (n_scenarios, n_locations, n_years) = size(final_arr)

    axlist::Tuple = (
        Dim{:timestep}(range(first_year, first_year + n_years - 1)),
        Dim{:location}(range(1, n_locations)),
        Dim{:scenario}(range(1, n_scenarios)),
    )

    props::Dict{String,Any} = Dict(
        "variable" => variable_name,
        "time units" => "years"
    )

    return YAXArray(axlist, permutedims(final_arr, [3, 2, 1]), props)
end
function _extract_variable(
    variable_arr::Vector{<:Array{<:Real, 4}}, variable_name::String, first_year::Int
)::YAXArray
    final_arr = _join_variable_arr(variable_arr)
    (n_scenarios, n_locations, n_years, n_groups) = size(final_arr)

    axlist::Tuple = (
        Dim{:timestep}(range(first_year, first_year + n_years - 1)),
        Dim{:location}(range(1, n_locations)),
        Dim{:group}(range(1, n_groups)),
        Dim{:scenario}(range(1, n_scenarios)),
    )

    props::Dict{String,Any} = Dict(
        "variable" => variable_name,
        "time units" => "years"
    )

    return YAXArray(axlist, permutedims(final_arr, [3, 2, 4, 1]), props)
end

"""
    to_NetCDF(
        filepath::String, 
        out_path::String, 
        variable_keys::String, 
        variable_names::Vector{Symbol}
    )::Nothing
"""
function to_NetCDF(
    directory::String,
    climate_scenario::String,
    out_path::String,
    variable_keys::Vector{String},
    variable_names::Vector{Symbol}
)::Nothing
    filenames::Vector{String} = _get_matfiles(directory, climate_scenario)
    matfiles = [read(matopen(fn)) for fn in filenames]

    arr_list = Vector{YAXArray}(undef, length(variable_keys) + 4)
    first_year::Float64 = matfiles[1]["YEARS"][1, 1]
    @info "Creating YAXArrays"
    for (ind, key) in enumerate(variable_keys)
        arr_list[ind] = _extract_variable(matfiles, key, first_year)
    end

    arr_list[length(variable_keys) + 1] = _extract_nongrazable_area(matfiles)
    arr_list[length(variable_keys) + 2] = _extract_reef_area(matfiles[1]) # asuumed duplicates
    arr_list[length(variable_keys) + 3] = _extract_lat(matfiles[1])
    arr_list[length(variable_keys) + 4] = _extract_lon(matfiles[1])

    props = _construct_properties(matfiles, filenames, climate_scenario)

    final_dataset = YAXArrays.Dataset(; (variable_names .=> arr_list)..., properties=props)
    @info "Saving YAXArrays"
    savedataset(final_dataset, path=out_path, driver=:netcdf, overwrite=true, compress=7)

    return nothing
end

"""
    function to_NetCDF_seperate(
        directory::String, climate_scenario::String;
        input_write_path::String="", output_write_path::String="", overwrite=false
    )::Nothing

Extract data from the given matfile into two NetCDF files, one for input variable and another for output variables.
"""
function to_NetCDF_seperate(
    directory::String, climate_scenario::String;
    input_write_path::String="", output_write_path::String="", overwrite=false
)::Nothing

    if input_write_path == ""
        input_write_path = "ReefMod_RCP" * "-inputs.nc"
    end

    if output_write_path == ""
        output_write_path = "ReefMod_RCP" * "-outputs.nc"
    end

    input_keys::Vector{String} = Vector{String}([
        "record_applied_DHWs",
        "record_applied_cyclones",
        "record_applied_bleaching_mortality",
        "coral_HT_mean",
        "coral_HT_var",
        "COTS_juv_densities",
        "COTS_adult_densities"
    ])

    input_var_names::Vector{Symbol} = Vector{Symbol}([
        :record_applied_DHWs,
        :record_applied_cyclones,
        :record_applied_bleaching_mortality,
        :coral_HT_mean,
        :coral_HT_var,
        :COTS_juv_densities,
        :COTS_adult_densities,
        :nongrazable,
        :reef_area,
        :lat,
        :lon,
    ])

    output_keys::Vector{String} = Vector{String}([
        "coral_larval_supply",
        "coral_cover_per_taxa",
        "nb_coral_offspring",
        "nb_coral_adol",
        "COTS_larval_supply",
        "nb_coral_adult",
        "nb_coral_juv",
        "nb_coral_recruit",
        "reef_shelter_volume_relative"
    ])

    output_var_names::Vector{Symbol} = Vector{Symbol}([
        :coral_larval_supply,
        :coral_cover_per_taxa,
        :nb_coral_offspring,
        :nb_coral_adol,
        :COTS_larval_supply,
        :nb_coral_adult,
        :nb_coral_juv,
        :nb_coral_recruit,
        :reef_shelter_volume_relative,
        :nongrazable,
        :reef_area,
        :lat,
        :lon,
    ])
    to_NetCDF(directory, climate_scenario, input_write_path, input_keys, input_var_names)
    to_NetCDF(directory, climate_scenario, output_write_path, output_keys, output_var_names)
    return nothing
end

"""
    function to_NetCDF_combined(
        directory::String, climate_scenario::String; out_path::String=""
    )::Nothing

Extract data from the given matfile into a single netcdf file containing all variables.
"""
function to_NetCDF_combined(
    directory::String, climate_scenario::String; out_path::String=""
)::Nothing

    if out_path == ""
        out_path = "ReefMod_RCP" * climate_scenario * ".nc"
    end

    var_keys::Vector{String} = Vector{String}([
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
    ])

    key_list::Vector{Symbol} = Vector{Symbol}([
        :record_applied_DHWs,
        :record_applied_cyclones,
        :record_applied_bleaching_mortality,
        :coral_HT_mean,
        :coral_HT_var,
        :COTS_juv_densities,
        :COTS_adult_densities,
        :coral_larval_supply,
        :coral_cover_per_taxa,
        :nb_coral_offspring,
        :nb_coral_adol,
        :COTS_larval_supply,
        :nb_coral_adult,
        :nb_coral_juv,
        :nb_coral_recruit,
        :reef_shelter_volume_relative,
        :nongrazable,
        :reef_area,
        :lat,
        :lon,
    ])

    to_NetCDF(directory, climate_scenario, out_path, var_keys, key_list)
    return nothing
end
