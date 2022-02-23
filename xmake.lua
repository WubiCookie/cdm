add_rules("mode.release", "mode.debug", "mode.releasedbg")
add_requires("catch2")

set_arch("x64")
set_plat("windows")

-- Generates a hash key made of packages confs/version, for CI
task("dephash")
	on_run(function ()
		import("core.project.project")
		import("private.action.require.impl.package")

		local requires, requires_extra = project.requires_str()

		local key = {}
		for _, instance in irpairs(package.load_packages(requires, {requires_extra = requires_extra})) do
			table.insert(key, instance:name() .. "-" .. instance:version_str() .. "-" .. instance:buildhash())
		end

		table.sort(key)

		key = table.concat(key, ",")
		print(hash.uuid4(key):gsub('-', ''):lower())
	end)

	set_menu {
		usage = "xmake dephash",
		description = "Outputs a hash key of current dependencies version/configuration"
	}
task_end()

target("testCdm")
	set_kind("binary")
	set_languages("c++20")
	add_files("tests/*.cpp")
	add_headerfiles("*.h*")
	add_packages("catch2")
	add_includedirs(".", "tests")
