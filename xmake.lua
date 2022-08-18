add_rules("mode.release", "mode.debug", "mode.releasedbg")
add_requires("catch2 2.13.9")

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

local tests = {
	"direction",
	"line",
	"matrix3",
	"matrix4",
	"misc",
	"perspective",
	"plane",
	"quaternion",
	"segment3",
	"static_pi_fraction_cos",
	"static_pi_fraction_sin",
	"static_pi_fraction_tan",
	"transform3",
	"transform3__multiplication",
	-- "transform3__rotation",
	"transform3__translation",
	"value_domain",
	"vector2",
	"vector3",
	"vector4",
}

for _,v in pairs(tests) do
	target("cdm_test_"..v)
		set_kind("binary")
		set_languages("c++20")
		add_files("tests/"..v..".cpp")
		add_headerfiles("*.hpp")
		add_headerfiles("tests/*.hpp")
		add_packages("catch2")
		add_includedirs(".", "tests")
	target_end()
end
