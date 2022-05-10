# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT license.

cmake_minimum_required(VERSION 3.12)

project(TestVectors VERSION 3.6.6 LANGUAGES CXX)

if(NOT DEFINED SEAL_BUILD_TESTVECTORS)
    message("yes")
    set(SEAL_BUILD_TESTVECTORS ON)

    # Import Microsoft SEAL
    find_package(SEAL 3.6.6 EXACT REQUIRED)

    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
endif()
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)

if(SEAL_BUILD_TESTVECTORS)
    message("yes again")
    add_executable(testvectors)

    target_sources(testvectors
        PRIVATE
            ${CMAKE_CURRENT_LIST_DIR}/testvectors.cpp
    )

    if(TARGET SEAL::seal)
        target_link_libraries(testvectors PRIVATE SEAL::seal)
    elseif(TARGET SEAL::seal_shared)
        target_link_libraries(testvectors PRIVATE SEAL::seal_shared)
    else()
        message(FATAL_ERROR "Cannot find target SEAL::seal or SEAL::seal_shared")
    endif()
endif()
