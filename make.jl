using PackageCompiler
if length(ARGS)>0 && ARGS[1] == "precompile"
    create_app(
        ".", 
        "compiled1"; 
        precompile_execution_file=realpath("precompile.jl"), 
        executables=["fastaExt" => "fastaext", 
                    "build_index" => "build_index", 
                    "build_header_index" => "build_header_index", 
                    "extgen" => "extgen"], 
        force=true, 
        cpu_target="native", 
        filter_stdlibs=true)
else
    create_app(
        ".", 
        "compiled"; 
        executables=["fastaExt" => "fastaext", 
                    "build_index" => "build_index", 
                    "build_header_index" => "build_header_index", 
                    "extgen" => "extgen"], 
        force=true, 
        cpu_target="native", 
        filter_stdlibs=true)
end
