
add_rules("mode.debug", "mode.release")
add_requires("eigen")

set_policy("build.warning", true)
set_warnings("all", "extra")

-- -- -- -- -- -- -- -- -- -- 

target("geok")
    set_kind("binary")
    set_languages("c++20")
    set_optimize("fastest")
    add_includedirs("./src", "./src/cddlib")
    add_deps("geoklib")
    add_files("./tests/*.cpp")
    add_packages("eigen")

target("geoklib")
    set_kind("static")
    set_languages("c11")
    set_optimize("fastest")
    add_includedirs("./src/cddlib")
    add_files("./src/cddlib/*.c")

