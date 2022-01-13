add_rules("mode.debug", "mode.release")
add_requires("catch2")
set_languages("c++17")

target("testCdm")
    set_kind("binary")
    add_files("main.cpp", "tests.cpp")
    add_packages("catch2")
    add_includedirs("..")