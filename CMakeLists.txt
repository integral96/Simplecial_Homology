cmake_minimum_required(VERSION 3.5)

project(my_yemp ${PROJECT_NAME}  LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

ADD_DEFINITIONS(-DQT_NO_KEYWORDS)
find_package (Boost ${BOOST_MIN_VERSION} COMPONENTS context
    filesystem date_time timer fiber
    thread system program_options regex log log_setup)

find_package(TBB COMPONENTS tbb)



set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -ltbb -lboost_log_setup -pthread --pedantic")
#set(CMAKE_CXX_FLAGS "-Wall -Wextra -fPIC -DBOOST_LOG_DYN_LINK")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftemplate-depth=2800")

#if(CMAKE_EXPORT_COMPILE_COMMANDS)
#    set(CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES ${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES})
#endif()

file(GLOB SOURCES
    *.hpp
    *.h
    *.cpp
)

add_executable(${PROJECT_NAME} ${SOURCES} )
target_link_libraries (${PROJECT_NAME} TBB::tbb ${Boost_LIBRARIES}

    )

