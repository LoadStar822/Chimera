if(NOT TARGET cereal::cereal)
    add_library(cereal::cereal INTERFACE IMPORTED)

    get_filename_component(_cereal_config_dir "${CMAKE_CURRENT_LIST_FILE}" PATH)
    get_filename_component(_cereal_project_root "${_cereal_config_dir}/../.." ABSOLUTE)
    set(_cereal_include_dir "${_cereal_project_root}/include/seqan3/submodules/cereal/include")

    set_target_properties(cereal::cereal PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${_cereal_include_dir}"
    )
endif()
