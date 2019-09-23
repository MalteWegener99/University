solution "Eigen Test"
    language "C++"
    buildoptions "-std=c++17"
    
    configurations {"Debug", "Release", "Release.Arm"}

    configuration { "configurations:Debug"}
        flags {"Symbols"}

    configuration {"configurations:Release"}
        defines {"NDEBUG"}
        flags {"Optimize"}
    configuration {"configurations:Release.Arm"}
    	architecture {"ARM"}
	defines {"NDEBUG"}
	flags {"Optimize"}

    targetdir ("Build/Bin")
    objdir ("Build/Obj")

    project "App"
        kind "ConsoleApp"
        language "C++"
        files {"Projects/App/**"}

        includedirs {"../../../../../../../../usr/include/eigen3", "../../../../../../../../usr/include/eigen3/unsupported"}

