set_project("cdm")

add_rules("mode.release", "mode.debug", "mode.releasedbg")

set_arch("x64")

target("cdm_maths_tests")
	-- set_default(false)
	set_kind("binary")
	set_languages("c++17")
	add_files("tests/cdm_maths_tests.cpp")
	add_headerfiles("*.h*")
	add_includedirs(".")
