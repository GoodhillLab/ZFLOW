CMAKE_MINIMUM_REQUIRED(VERSION 3.15.0)
PROJECT(zebrafish3d)

# List your source files here - this example has just one.
SET(SOURCE_FILES example.cpp IBZebraFishKinematics3d.cpp)
ADD_EXECUTABLE(main3d ${SOURCE_FILES})

FIND_PACKAGE(IBAMR 0.14.0 REQUIRED)
TARGET_LINK_LIBRARIES(main3d IBAMR::IBAMR3d)
# IBAMR saves the flags it used to compile - you can reuse them if you want to
SET(CMAKE_CXX_FLAGS ${IBAMR_CXX_FLAGS})
