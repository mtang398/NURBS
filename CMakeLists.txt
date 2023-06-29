cmake_minimum_required(VERSION 3.14)  # FetchContent is available in 3.14+
project(NURBS)

include(FetchContent)

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_VIEWER            "Use OpenGL viewer"  ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"          ON)

# Provide hint for finding BLAS
set(BLA_VENDOR OpenBLAS)
find_package(BLAS REQUIRED)

FetchContent_Declare(
  libigl
  GIT_REPOSITORY https://github.com/libigl/libigl.git
  GIT_TAG        v2.3.0  # adjust with the version you want
)

FetchContent_GetProperties(libigl)
if(NOT libigl_POPULATED)
  FetchContent_Populate(libigl)
  add_subdirectory(${libigl_SOURCE_DIR} ${libigl_BINARY_DIR})
endif()

# Add your project files
file(GLOB SRCFILES *.cpp)
add_executable(${PROJECT_NAME}_bin ${SRCFILES})
target_link_libraries(${PROJECT_NAME}_bin igl::core igl::opengl_glfw igl::opengl_glfw_imgui)
