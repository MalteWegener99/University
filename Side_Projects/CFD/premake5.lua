workspace "CFD"
   configurations { "Debug", "Release" }

project "Main"
   kind "ConsoleApp"
   language "C++"
   cppdialect "c++17"

   targetdir "bin/%{cfg.buildcfg}"
   includedirs {
       "../../../../../../../../../usr/include/eigen3",
       "../../../../../../../../../usr/include/eigen3/unsupported"
   }

   files { "**.h", "**.c" }

   filter "configurations:Debug"
      defines { "DEBUG" }
      symbols "On"

   filter "configurations:Release"
      defines { "NDEBUG" }
      optimize "On"