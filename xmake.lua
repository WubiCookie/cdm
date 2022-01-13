add_rules("mode.release", "mode.debug", "mode.releasedbg")
add_requires("catch2")

set_arch("x64")

target("testCdm")
	set_kind("binary")
	set_languages("c++17")
	add_files("tests/*.cpp")
	add_headerfiles("*.h*")
	add_packages("catch2")
	add_includedirs(".")
